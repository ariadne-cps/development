/***************************************************************************
 *            concurrency/task_runner.tpl.hpp
 *
 *  Copyright  2007-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file concurrency/task_runner.tpl.hpp
 *  \brief Implementation code for runner classes.
 */

#ifndef ARIADNE_TASK_RUNNER_TPL_HPP
#define ARIADNE_TASK_RUNNER_TPL_HPP

#include "task.tpl.hpp"
#include "task_runner.hpp"
#include "task_interface.hpp"
#include "concurrency_manager.hpp"
#include "configuration/configurable.tpl.hpp"

namespace Ariadne {

template<class C> TaskRunnable<C>::TaskRunnable(ConfigurationType const& configuration) : Configurable<C>(configuration) {
    ConcurrencyManager::instance().choose_runner_for(*this);
}

template<class C> void TaskRunnable<C>::set_runner(SharedPointer<TaskRunnerInterface<C>> const& runner) {
    this->_runner = runner;
}

template<class C> SharedPointer<TaskRunnerInterface<C>>& TaskRunnable<C>::runner() {
    return _runner;
}

template<class C> SharedPointer<TaskRunnerInterface<C>> const& TaskRunnable<C>::runner() const {
    return _runner;
}

template<class C> class TaskRunnerBase : public TaskRunnerInterface<C> {
  public:
    typedef typename TaskRunnerInterface<C>::TaskType TaskType;
    typedef typename TaskRunnerInterface<C>::InputType InputType;
    typedef typename TaskRunnerInterface<C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<C>::ConfigurationType ConfigurationType;

    TaskRunnerBase(ConfigurationType const& configuration)
        : _task(new TaskType()), _configuration(configuration) { }

    TaskType& task() override { return *_task; };
    TaskType const& task() const override { return *_task; };
    ConfigurationType const& configuration() const override { return _configuration; }

    virtual ~TaskRunnerBase() = default;

  protected:
    SharedPointer<TaskType> const _task;
    ConfigurationType const _configuration;
};

template<class C> SequentialRunner<C>::SequentialRunner(ConfigurationType const& configuration) : TaskRunnerBase<C>(configuration) { }

template<class C> void SequentialRunner<C>::push(InputType const& input) {
    OutputType result = this->_task->run(input,this->configuration());
    auto failed_constraints = this->_task->ranking_space().failed_critical_constraints(input,result);
    if (not failed_constraints.empty()) throw CriticalRankingFailureException<C>(failed_constraints);
    _last_output.reset(new OutputType(result));
}

template<class C> auto SequentialRunner<C>::pull() -> OutputType {
    return *_last_output;
}

template<class C> void DetachedRunner<C>::_loop() {
    while(true) {
        std::unique_lock<std::mutex> locker(_input_mutex);
        _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
        if (_terminate) break;
        auto input = _input_buffer.pop();
        OutputType output = this->_task->run(input,this->configuration());
        _output_buffer.push(output);
        _output_availability.notify_all();
    }
}

template<class C> DetachedRunner<C>::DetachedRunner(ConfigurationType const& configuration)
        : TaskRunnerBase<C>(configuration),
          _thread(this->_task->name(), [this]() { _loop(); }),
          _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
          _last_used_input(1), _active(false), _terminate(false) { }

template<class C> DetachedRunner<C>::~DetachedRunner() {
    _terminate = true;
    _input_availability.notify_all();
}

template<class C> void DetachedRunner<C>::push(InputType const& input) {
    if (not _active) {
        _active = true;
        _thread.activate();
    }
    _input_buffer.push(input);
    _last_used_input.push(input);
    _input_availability.notify_all();
}

template<class C> auto DetachedRunner<C>::pull() -> OutputType {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>0; });
    auto result = _output_buffer.pop();
    auto failed_constraints = this->_task->ranking_space().failed_critical_constraints(_last_used_input.pop(),result);
    if (not failed_constraints.empty()) throw CriticalRankingFailureException<C>(failed_constraints);
    return result;
}

template<class O> class ParameterSearchOutputBufferData {
  public:
    ParameterSearchOutputBufferData(O const& output, DurationType const& execution_time, ConfigurationSearchPoint const& point) : _output(output), _execution_time(execution_time), _point(point) { }
    ParameterSearchOutputBufferData(ParameterSearchOutputBufferData<O> const& p) : _output(p._output), _execution_time(p._execution_time), _point(p._point) { }
    ParameterSearchOutputBufferData& operator=(ParameterSearchOutputBufferData<O> const& p) {
        _output = p._output;
        _execution_time = p._execution_time;
        _point = p._point;
        return *this;
    };
    O const& output() const { return _output; }
    DurationType const& execution_time() const { return _execution_time; }
    ConfigurationSearchPoint const& point() const { return _point; }
  private:
    O _output;
    DurationType _execution_time;
    ConfigurationSearchPoint _point;
};

template<class C> void ParameterSearchRunner<C>::_loop() {
    while(true) {
        std::unique_lock<std::mutex> locker(_input_mutex);
        _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
        locker.unlock();
        if (_terminate) break;
        if (_input_buffer.size() == 0) std::cout << "Input buffer is empty." << std::endl;
        auto pkg = _input_buffer.pop();
        auto cfg = make_singleton(this->configuration(),pkg.second);
        try {
            auto start = std::chrono::high_resolution_clock::now();
            auto output = this->_task->run(pkg.first,cfg);
            auto end = std::chrono::high_resolution_clock::now();
            auto execution_time = std::chrono::duration_cast<DurationType>(end-start);
            ARIADNE_LOG_PRINTLN("task for " << pkg.second << " completed in " << execution_time.count() << " us");
            _output_buffer.push(OutputBufferContentType(output,execution_time,pkg.second));
        } catch (std::exception& e) {
            ++_failures;
            ARIADNE_LOG_PRINTLN("task failed: " << e.what());
        }
        _output_availability.notify_all();
    }
}

template<class C> ParameterSearchRunner<C>::ParameterSearchRunner(ConfigurationType const& configuration, Nat concurrency)
        : TaskRunnerBase<C>(configuration), _concurrency(concurrency),
          _failures(0), _last_used_input(1),
          _input_buffer(InputBufferType(concurrency)),_output_buffer(OutputBufferType(concurrency)),
          _active(false), _terminate(false) {
    for (Nat i=0; i<concurrency; ++i)
        _threads.append(SharedPointer<LoggableSmartThread>(new LoggableSmartThread(this->_task->name() + (concurrency>=10 and i<10 ? "0" : "") + to_string(i), [this]() { _loop(); })));
}

template<class C> ParameterSearchRunner<C>::~ParameterSearchRunner() {
    _terminate = true;
    _input_availability.notify_all();
}

template<class C> void ParameterSearchRunner<C>::push(InputType const& input) {
    if (not _active) {
        _active = true;
        auto shifted = this->_configuration.search_space().initial_point().make_random_shifted(_concurrency);
        for (auto point : shifted) _points.push(point);
        for (auto thread : _threads) thread->activate();
    }
    for (SizeType i=0; i<_concurrency; ++i) {
        _input_buffer.push({input,_points.front()});
        _points.pop();
    }
    _last_used_input.push(input);
    _input_availability.notify_all();
}

template<class C> auto ParameterSearchRunner<C>::pull() -> OutputType {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>=_concurrency-_failures; });
    ARIADNE_LOG_PRINTLN("received " << _concurrency-_failures << " completed tasks");
    _failures=0;

    InputType input = _last_used_input.pop();
    Map<ConfigurationSearchPoint,Pair<OutputType,DurationType>> outputs;
    while (_output_buffer.size() > 0) {
        auto io_data = _output_buffer.pop();
        outputs.insert(Pair<ConfigurationSearchPoint,Pair<OutputType,DurationType>>(
                io_data.point(),{io_data.output(),io_data.execution_time()}));
    }
    auto rankings = this->_task->rank(outputs,input);
    ARIADNE_LOG_PRINTLN_VAR(rankings);

    Set<ConfigurationSearchPoint> new_points;
    SizeType cnt = 0;
    for (auto it = rankings.rbegin(); it != rankings.rend(); ++it) {
        new_points.insert(it->point());
        ++cnt;
        if (cnt >= _concurrency/2) break;
    }
    new_points = make_extended_set_by_shifting(new_points, _concurrency);
    for (auto p : new_points) _points.push(p);
    ARIADNE_LOG_PRINTLN_VAR(new_points);

    auto best = rankings.rbegin()->point();
    if (rankings.rbegin()->critical_failures() > 0) {
        throw CriticalRankingFailureException<C>(this->_task->ranking_space().failed_critical_constraints(input, outputs.get(best).first));
    }
    ConcurrencyManager::instance().append_best_ranking(*rankings.rbegin());
    auto best_output = outputs.get(best).first;

    return best_output;
}

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_TPL_HPP
