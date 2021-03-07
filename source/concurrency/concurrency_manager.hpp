/***************************************************************************
 *            concurrency/concurrency_manager.hpp
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

/*! \file concurrency/concurrency_manager.hpp
 *  \brief Singleton class for managing concurrency across the library
 */

#ifndef ARIADNE_CONCURRENCY_MANAGER_HPP
#define ARIADNE_CONCURRENCY_MANAGER_HPP

#include <thread>
#include <algorithm>
#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../concurrency/task_runner.hpp"
#include "../concurrency/task_execution_ranking.hpp"
#include "../configuration/configuration_property_interface.hpp"

namespace Ariadne {

//! \brief Manages threads and sets runners based on concurrency availability.
class ConcurrencyManager {
    typedef Map<ConfigurationPropertyPath,List<ExactDouble>> PropertyRefinementsMap;
  private:
    ConcurrencyManager();
  public:
    ConcurrencyManager(ConcurrencyManager const&) = delete;
    void operator=(ConcurrencyManager const&) = delete;

    static ConcurrencyManager& instance() {
        static ConcurrencyManager instance;
        return instance;
    }

    //! \brief Choose the proper runner for \a runnable
    template<class T> void choose_runner_for(TaskRunnable<T>& runnable) const {
        SharedPointer<TaskRunnerInterface<T>> runner;
        auto const& cfg = runnable.configuration();
        if (_concurrency > 1 and not cfg.is_singleton())
            runner.reset(new ParameterSearchRunner<T>(cfg,std::min(_concurrency,cfg.search_space().total_points())));
        else if (_concurrency == 1 and not cfg.is_singleton()) {
            auto point = cfg.search_space().initial_point();
            ARIADNE_LOG_PRINTLN_AT(1,"The configuration is not singleton: using point " << point << " for sequential running.");
            runner.reset(new SequentialRunner<T>(make_singleton(cfg,point)));
        } else
            runner.reset(new SequentialRunner<T>(cfg));
        runnable.set_runner(runner);
    }

    SizeType maximum_concurrency() const;
    SizeType concurrency() const;

    void set_concurrency(SizeType value);

    //! \brief The best points saved
    List<TaskExecutionRanking> best_rankings() const;
    void append_best_ranking(TaskExecutionRanking const& point);
    void clear_best_rankings();

    //! \brief Print best rankings in a .m file for plotting
    void print_best_rankings() const;

    //! \brief Return the optimal point (i.e., the most common value for all dimensions)
    List<int> optimal_point() const;

  private:
    const SizeType _maximum_concurrency;
    SizeType _concurrency;
    std::mutex _data_mutex;
    List<TaskExecutionRanking> _best_rankings;
};

} // namespace Ariadne

#endif // ARIADNE_CONCURRENCY_MANAGER_HPP
