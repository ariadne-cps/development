/***************************************************************************
 *            concurrency/task_runner.hpp
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

/*! \file concurrency/task_runner.hpp
 *  \brief Runner classes.
 */

#ifndef ARIADNE_TASK_RUNNER_HPP
#define ARIADNE_TASK_RUNNER_HPP

#include "loggable_smart_thread.hpp"
#include "buffer.hpp"
#include "task_runner_interface.hpp"
#include "configuration/configuration_search_point.hpp"
#include "configuration/configuration_search_space.hpp"
#include "task_execution_ranking.hpp"
#include "configuration/configuration_property.tpl.hpp"

namespace Ariadne {

template<class C> class TaskRunnerBase;

//! \brief Run a task sequentially.
//! \details Used to provide a sequential alternative to any thread-based implementation.
template<class C>
class SequentialRunner final : public TaskRunnerBase<C> {
    friend class ConcurrencyManager;
    typedef typename TaskRunnerBase<C>::InputType InputType;
    typedef typename TaskRunnerBase<C>::OutputType OutputType;
    typedef typename TaskRunnerBase<C>::ConfigurationType ConfigurationType;
  protected:
    SequentialRunner(ConfigurationType const& configuration);
  public:
    Void push(InputType const& input) override final;
    OutputType pull() override final;

private:
    SharedPointer<OutputType> _last_output;
};

//! \brief Run a task in a detached thread, allowing other processing between pushing and pulling.
template<class C>
class DetachedRunner final : public TaskRunnerBase<C> {
    friend class ConcurrencyManager;
    typedef typename TaskRunnerBase<C>::InputType InputType;
    typedef typename TaskRunnerBase<C>::OutputType OutputType;
    typedef typename TaskRunnerBase<C>::ConfigurationType ConfigurationType;
    typedef Buffer<InputType> InputBufferType;
    typedef Buffer<OutputType> OutputBufferType;
  protected:
    DetachedRunner(ConfigurationType const& configuration);
  public:
    virtual ~DetachedRunner();

    Void push(InputType const& input) override final;
    OutputType pull() override final;

private:
    Void _loop();
private:
    LoggableSmartThread _thread;
    InputBufferType _input_buffer;
    OutputBufferType _output_buffer;
    Buffer<InputType> _last_used_input;
    std::atomic<bool> _active;
    std::atomic<bool> _terminate;
    std::mutex _input_mutex;
    std::condition_variable _input_availability;
    std::mutex _output_mutex;
    std::condition_variable _output_availability;
};

typedef std::chrono::microseconds DurationType;

template<class O> class ParameterSearchOutputBufferData;

//! \brief Run a task by detached concurrent search into the parameter space.
template<class C> class ParameterSearchRunner final : public TaskRunnerBase<C> {
    friend class ConcurrencyManager;
    typedef typename TaskRunnerBase<C>::InputType InputType;
    typedef typename TaskRunnerBase<C>::OutputType OutputType;
    typedef typename TaskRunnerBase<C>::ConfigurationType ConfigurationType;
    typedef Pair<InputType,ConfigurationSearchPoint> InputBufferContentType;
    typedef ParameterSearchOutputBufferData<OutputType> OutputBufferContentType;
    typedef Buffer<InputBufferContentType> InputBufferType;
    typedef Buffer<OutputBufferContentType> OutputBufferType;
  protected:
    ParameterSearchRunner(ConfigurationType const& configuration, Nat concurrency);
  public:
    virtual ~ParameterSearchRunner();

    Void push(InputType const& input) override final;
    OutputType pull() override final;

private:
    Void _loop();
private:
    Nat const _concurrency; // Number of threads to be used
    std::atomic<Nat> _failures; // Number of failures after a given push, reset during pulling
    Buffer<InputType> _last_used_input;
    std::queue<ConfigurationSearchPoint> _points;
    // Synchronization
    List<SharedPointer<LoggableSmartThread>> _threads;
    InputBufferType _input_buffer;
    OutputBufferType _output_buffer;
    std::atomic<bool> _active;
    std::atomic<bool> _terminate;
    std::mutex _input_mutex;
    std::condition_variable _input_availability;
    std::mutex _output_mutex;
    std::condition_variable _output_availability;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_HPP
