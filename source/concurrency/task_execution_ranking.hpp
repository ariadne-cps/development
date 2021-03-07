/***************************************************************************
 *            concurrency/task_execution_ranking.hpp
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

/*! \file concurrency/task_execution_ranking.hpp
 *  \brief Class for defining the score for a point.
 */

#ifndef ARIADNE_TASK_EXECUTION_RANKING_HPP
#define ARIADNE_TASK_EXECUTION_RANKING_HPP

#include "configuration/configuration_search_point.hpp"
#include "../concurrency/task_ranking_parameter.hpp"

namespace Ariadne {

typedef double ScoreType;

class TaskExecutionRanking;

template<class R> class CriticalRankingFailureException : public std::runtime_error {
public:
    CriticalRankingFailureException(List<TaskRankingParameter<R>> const& parameters) : std::runtime_error("The execution has critical failures for these parameters: " + to_string(parameters)) { }
};

class TaskExecutionRanking : public WritableInterface {
public:
    TaskExecutionRanking(ConfigurationSearchPoint const& p, ScoreType const& s, SizeType const& permissive_failures, SizeType const& critical_failures);
    ConfigurationSearchPoint const& point() const;
    ScoreType const& score() const;
    SizeType const& permissive_failures() const;
    SizeType const& critical_failures() const;
    //! \brief Ordering is based on number of failures, followed by cost
    Bool operator<(TaskExecutionRanking const& s) const;

    virtual OutputStream& _write(OutputStream& os) const;
private:
    ConfigurationSearchPoint _point;
    ScoreType _score;
    SizeType _permissive_failures;
    SizeType _critical_failures;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_EXECUTION_RANKING_HPP
