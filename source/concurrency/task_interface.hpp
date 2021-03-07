/***************************************************************************
 *            concurrency/task_interface.hpp
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

/*! \file concurrency/task_interface.hpp
 *  \brief The interface for tasks.
 */

#ifndef ARIADNE_TASK_INTERFACE_HPP
#define ARIADNE_TASK_INTERFACE_HPP

#include <chrono>
#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../utility/string.hpp"

namespace Ariadne {

class ConfigurationSearchPoint;
class TaskExecutionRanking;
class ConfigurationSearchSpace;
template<class R> class TaskRankingSpace;

typedef std::chrono::microseconds DurationType;

template<class R> struct TaskInput;
template<class R> struct TaskOutput;
template<class R> struct TaskObjective;
template<class R> struct Task;
template<class R> struct Configuration;

template<class R>
class TaskInterface {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef Configuration<R> ConfigurationType;

    //! \brief The name of the task, to be used for thread naming
    virtual String name() const = 0;
    //! \brief Return the ranking space for the task
    virtual TaskRankingSpace<R> const& ranking_space() const = 0;
    //! \brief Set the ranking space for the task
    virtual Void set_ranking_space(TaskRankingSpace<R> const& space) = 0;

    //! \brief The task to be performed, taking \a in as input and \a cfg as a configuration of the parameters
    virtual OutputType run(InputType const& in, ConfigurationType const& cfg) const = 0;
    //! \brief Evaluate the costs of points from output and execution time, possibly using the input \a in
    virtual Set<TaskExecutionRanking> rank(Map<ConfigurationSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& in) const = 0;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_INTERFACE_HPP
