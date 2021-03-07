/***************************************************************************
 *            concurrency/task.tpl.hpp
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

/*! \file concurrency/task.tpl.hpp
 *  \brief Base code for a task.
 */

#ifndef ARIADNE_TASK_TPL_HPP
#define ARIADNE_TASK_TPL_HPP

#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../utility/string.hpp"
#include "task_interface.hpp"
#include "task_ranking_space.hpp"
#include "configuration/configuration_search_space.hpp"

namespace Ariadne {

class ConfigurationSearchPoint;
class TaskExecutionRanking;
class ConfigurationSearchSpace;

//! \brief The base for parameter search tasks
//! \details Useful to streamline task construction
template<class R>
class ParameterSearchTaskBase : public TaskInterface<R> {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
  protected:
    ParameterSearchTaskBase(String const& name)
        : _name(name), _ranking_space(SharedPointer<TaskRankingSpace<R>>(TaskRankingSpaceBuilder<R>().build().clone())) {}
  public:
    String name() const override { return _name; }
    TaskRankingSpace<R> const& ranking_space() const override { return *_ranking_space; }
    Void set_ranking_space(TaskRankingSpace<R> const& space) override { _ranking_space.reset(space.clone()); }

    Set<TaskExecutionRanking> rank(Map<ConfigurationSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& input) const override { return _ranking_space->rank(data, input); }

  private:
    String const _name;
    SharedPointer<TaskRankingSpace<R>> _ranking_space;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_TPL_HPP
