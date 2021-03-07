/***************************************************************************
 *            concurrency/task_execution_ranking.cpp
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

#include "../concurrency/task_execution_ranking.hpp"

namespace Ariadne {

TaskExecutionRanking::TaskExecutionRanking(ConfigurationSearchPoint const& p,
                                           ScoreType const& s,
                                           SizeType const& permissive_failures,
                                           SizeType const& critical_failures)
           : _point(p), _score(s), _permissive_failures(permissive_failures), _critical_failures(critical_failures) { }

ConfigurationSearchPoint const&
TaskExecutionRanking::point() const {
    return _point;
}

ScoreType const&
TaskExecutionRanking::score() const {
    return _score;
}

SizeType const&
TaskExecutionRanking::permissive_failures() const {
    return _permissive_failures;
}

SizeType const&
TaskExecutionRanking::critical_failures() const {
    return _critical_failures;
}

Bool TaskExecutionRanking::operator<(TaskExecutionRanking const& s) const {
    if (this->_critical_failures > s._critical_failures) return true;
    else if (this->_critical_failures < s._critical_failures) return false;
    else if (this->_permissive_failures > s._permissive_failures) return true;
    else if (this->_permissive_failures < s._permissive_failures) return false;
    else return this->_score < s._score;
}

OutputStream& TaskExecutionRanking::_write(OutputStream& os) const {
    return os << "{" << _point << ":" << _score << (_permissive_failures > 0 ? (",P:" + to_string(_permissive_failures)) : "")
              << (_critical_failures > 0 ? (",C:" + to_string(_critical_failures)) : "") << "}";
}
} // namespace Ariadne
