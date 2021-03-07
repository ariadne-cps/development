/***************************************************************************
 *            concurrency/task_raking_parameter.hpp
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

/*! \file concurrency/task_ranking_parameter.hpp
 *  \brief Classes for handling parameters for ranking the results of a task.
 */

#ifndef ARIADNE_TASK_RANKING_PARAMETER_HPP
#define ARIADNE_TASK_RANKING_PARAMETER_HPP

#include <functional>
#include <chrono>
#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/handle.hpp"

namespace Ariadne {

template<class R> struct TaskInput;
template<class R> struct TaskOutput;
template<class R> struct TaskObjective;

enum class OptimisationCriterion { MINIMISE, MAXIMISE };
inline std::ostream& operator<<(std::ostream& os, const OptimisationCriterion opt) {
    switch (opt) {
        case OptimisationCriterion::MAXIMISE: os << "MAXIMISE"; break;
        case OptimisationCriterion::MINIMISE: os << "MINIMISE"; break;
        default: ARIADNE_FAIL_MSG("Unhandled OptimisationCriterion value.");
    }
    return os;
}

//! \brief Enumeration for the severity of satisfying a constraint
//! \details NONE: there actually is no constraint
//!          PERMISSIVE: satisfying the constraint is only desired
//!          CRITICAL: satisfying the constraint is mandatory
enum class RankingConstraintSeverity { NONE, PERMISSIVE, CRITICAL };
inline std::ostream& operator<<(std::ostream& os, const RankingConstraintSeverity severity) {
    switch (severity) {
        case RankingConstraintSeverity::NONE: os << "NONE"; break;
        case RankingConstraintSeverity::PERMISSIVE: os << "PERMISSIVE"; break;
        case RankingConstraintSeverity::CRITICAL: os << "CRITICAL"; break;
        default: ARIADNE_FAIL_MSG("Unhandled RankingConstraintSeverity value.");
    }
    return os;
}

typedef double ScoreType;
typedef std::chrono::microseconds DurationType;

template<class R> class TaskRankingParameterInterface : public WritableInterface {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;

    virtual String const& name() const = 0;
    virtual OptimisationCriterion optimisation() const = 0;
    virtual RankingConstraintSeverity severity() const = 0;
    virtual Bool is_scalar() const = 0;
    virtual Bool uses_objective() const = 0;
    virtual Bool discard(InputType const& input) const = 0;
    virtual ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const = 0;
    virtual ScoreType threshold(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const = 0;
    virtual SizeType dimension(InputType const& input) const = 0;

    virtual TaskRankingParameterInterface* clone() const = 0;
    virtual ~TaskRankingParameterInterface() = default;

    OutputStream& _write(OutputStream& os) const override { return os << *this; }

    friend OutputStream& operator<<(OutputStream& os, TaskRankingParameterInterface<R> const& p) {
        os << "{'" << p.name() << "'," << p.optimisation() << "," <<
           (p.is_scalar() ? "SCALAR":"VECTOR") << "," <<
           (p.uses_objective() ? "OBJECTIVE":"NO_OBJECTIVE") << "}"; return os; }
};

template<class R> class TaskRankingParameterBase : public TaskRankingParameterInterface<R> {
  public:
    TaskRankingParameterBase(String const& name, OptimisationCriterion const& opt, RankingConstraintSeverity const& severity)
            : _name(name), _optimisation(opt), _severity(severity) { }

    String const& name() const override { return _name; }
    OptimisationCriterion optimisation() const override { return _optimisation; }
    RankingConstraintSeverity severity() const override { return _severity; }

  private:
    String const _name;
    OptimisationCriterion const _optimisation;
    RankingConstraintSeverity const _severity;
};

template<class R> class ScalarRankingParameter : public TaskRankingParameterBase<R> {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    ScalarRankingParameter(String const& name, OptimisationCriterion const& opt, std::function<ScoreType(InputType const&, OutputType const&, DurationType const&)> const rfunc)
        : TaskRankingParameterBase<R>(name, opt, RankingConstraintSeverity::NONE), _rfunc(rfunc) { }

    Bool is_scalar() const override { return true; }
    Bool uses_objective() const override { return false; }
    Bool discard(InputType const& input) const override { return false; }
    SizeType dimension(InputType const& input) const override { return 1; }
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const override { return _rfunc(input, output, duration); }
    ScoreType threshold(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const override { ARIADNE_ERROR("Cannot compute threshold for non-objective scalar parameter"); return 0; }

    ScalarRankingParameter* clone() const override { return new ScalarRankingParameter(*this); }
  private:
    std::function<ScoreType(InputType const&, OutputType const&, DurationType const&)> const _rfunc;
};

template<class R> class ScalarObjectiveRankingParameter : public TaskRankingParameterBase<R> {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef TaskObjective<R> ObjectiveType;
    typedef std::function<ScoreType(InputType const&, OutputType const&, DurationType const&, ObjectiveType const&)> MeasureFunctionType;
    typedef std::function<Bool(InputType const&, ObjectiveType const&)> DiscardFunctionType;
    ScalarObjectiveRankingParameter(String const& name, OptimisationCriterion const& opt, RankingConstraintSeverity const& severity, ObjectiveType const& objective,
                                    MeasureFunctionType const& score, MeasureFunctionType const& threshold, DiscardFunctionType const& discard)
            : TaskRankingParameterBase<R>(name, opt, severity), _objective(objective),  _sfunc(score), _tfunc(threshold), _dfunc(discard) { }

    Bool is_scalar() const override { return true; }
    Bool uses_objective() const override { return true; }
    SizeType dimension(InputType const& input) const override { return 1; }
    Bool discard(InputType const& input) const override { return _dfunc(input,_objective); }
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const override {
        auto score = _sfunc(input, output, duration, _objective);
        return score;
    }
    ScoreType threshold(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const override { return _tfunc(input,output,duration,_objective); }

    ScalarObjectiveRankingParameter* clone() const override { return new ScalarObjectiveRankingParameter(*this); }
private:
    ObjectiveType const _objective;
    MeasureFunctionType const _sfunc;
    MeasureFunctionType const _tfunc;
    DiscardFunctionType const _dfunc;
};

template<class R> class VectorRankingParameter : public TaskRankingParameterBase<R> {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    VectorRankingParameter(String const& name, OptimisationCriterion const& opt, std::function<ScoreType(InputType const&, OutputType const&, DurationType const&, SizeType const&)> const rfunc, std::function<SizeType(InputType const&)> const dfunc)
            : TaskRankingParameterBase<R>(name, opt, RankingConstraintSeverity::NONE), _rfunc(rfunc), _dfunc(dfunc) { }

    Bool is_scalar() const override { return false; };
    Bool uses_objective() const override { return false; }

    SizeType dimension(InputType const& input) const override { return _dfunc(input); }
    Bool discard(InputType const& input) const override { return false; }
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx) const override { return _rfunc(input, output, duration, idx); }
    ScoreType threshold(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const override { ARIADNE_ERROR("Cannot compute threshold for non-objective vector parameter"); return 0; }

    VectorRankingParameter* clone() const override { return new VectorRankingParameter(*this); }
private:
    std::function<ScoreType(InputType const&, OutputType const&, DurationType const&, SizeType const&)> const _rfunc;
    std::function<SizeType(InputType const&)> const _dfunc;
};

template<class R> class TaskRankingParameter : public Handle<TaskRankingParameterInterface<R>> {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    using Handle<TaskRankingParameterInterface<R>>::Handle;
public:

    String const& name() const { return this->_ptr->name(); }
    OptimisationCriterion optimisation() const { return this->_ptr->optimisation(); }
    RankingConstraintSeverity severity() const { return this->_ptr->severity(); }
    Bool is_scalar() const { return this->_ptr->is_scalar(); };
    Bool uses_objective() const { return this->_ptr->uses_objective(); }
    Bool discard(InputType const& input) const {
        return this->_ptr->discard(input); }
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const {
        return this->_ptr->rank(input, output, duration, idx); }
    ScoreType threshold(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const {
        return this->_ptr->threshold(input, output, duration, idx); }

    SizeType dimension(InputType const& input) const { return this->_ptr->dimension(input); }
    friend OutputStream& operator<<(OutputStream& os, const TaskRankingParameter<R>& p) { return os << *p._ptr; }
};

//! \brief Template instance of the commonly-used execution time parameter for appraisal
template<class R> ScalarRankingParameter<R> execution_time_ranking("execution_time", OptimisationCriterion::MINIMISE, [](TaskInput<R> const& i, TaskOutput<R> const& o, DurationType const& d) { return d.count(); });

} // namespace Ariadne

#endif // ARIADNE_TASK_RANKING_PARAMETER_HPP
