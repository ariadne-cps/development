/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
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

#include "ariadne.hpp"
#include "output/progress_indicator.hpp"
#include "output/drawer.hpp"

using namespace Ariadne;

void print_all_rankings(List<List<TaskExecutionRanking>> const& rankings) {
    auto space = rankings.front().front().point().space();
    auto dimension = space.dimension();

    for (SizeType i = 0; i < dimension; ++i) {
        auto name = space.parameters()[i].path().last();
        std::ofstream file;
        file.open(name + ".m");
        file << "figure(1);\n";
        file << "xlim([0 630]);\n";
        file << "yticks(-100:1:100);\n";
        file << "hold on;\n";
        for (SizeType j = 0; j < rankings.size(); ++j) {
            auto size = rankings.at(j).size();
            file << "x = [1:" << size << "];\n";
            List<int> values;
            for (auto ranking : rankings.at(j)) values.push_back(ranking.point().coordinates()[i]);
            file << "y = " << values << ";\n";
            file << "plot(x,y,'ko');\n";
        }
        file << "hold off;\n";
        file.close();
    }

    std::ofstream file;
    file.open("soft_failures.m");
    file << "figure(1);\n";
    file << "xlim([0 630]);\n";
    file << "yticks(-100:1:100);\n";
    file << "hold on;\n";
    for (SizeType j = 0; j < rankings.size(); ++j) {
        auto size = rankings.at(j).size();
        file << "x = [1:" << size << "];\n";
        List<int> values;
        for (auto ranking : rankings.at(j)) values.push_back(ranking.permissive_failures());
        file << "y = " << values << ";\n";
        file << "plot(x,y,'ko');\n";
    }
    file << "hold off;\n";
    file.close();
}

int main(int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));
    Logger::instance().configuration().set_theme(TT_THEME_DARK);
    Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);
    ConcurrencyManager::instance().set_concurrency(5);
    Logger::instance().use_blocking_scheduler();

    ARIADNE_LOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField system({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    double max_err = 1e-8;
    auto sweeper1 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/10);
    auto sweeper2 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/100);
    auto sweeper3 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/1000);

    TaylorPicardIntegrator integrator(Configuration<TaylorPicardIntegrator>()
                                          .set_step_maximum_error(1e-8,1e-5)
                                          .set_maximum_temporal_order(9,12)
                                          .set_lipschitz_tolerance(0.05,0.5)
                                          .set_sweeper({sweeper1,sweeper2,sweeper3})
                                          );

    typedef VectorFieldEvolver E; typedef TaskInput<E> I; typedef TaskOutput<E> O; typedef TaskObjective<E> OBJ;

    E evolver(system,Configuration<E>().set_integrator(integrator));
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration());
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration().search_space());

    OBJ y_p275(y,PositiveFloatDPUpperBound(FloatDP(cast_exact(0.07),DoublePrecision())),Dyadic(cast_exact(6.48)));
    OBJ y_m275(y,PositiveFloatDPUpperBound(FloatDP(cast_exact(0.075),DoublePrecision())),Dyadic(cast_exact(3.15)));
    auto verification_p275 = ScalarObjectiveRankingParameter<E>(y.name(), OptimisationCriterion::MINIMISE, RankingConstraintSeverity::CRITICAL, y_p275,
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return o.reach.bounding_box()[obj.variable].upper_bound().get_d(); },
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return 2.75; },
                                      [](I const& i, OBJ const& obj) { return false; }
                                      );
    auto verification_m275 = ScalarObjectiveRankingParameter<E>(y.name(), OptimisationCriterion::MAXIMISE, RankingConstraintSeverity::CRITICAL, y_m275,
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return o.reach.bounding_box()[obj.variable].lower_bound().get_d(); },
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return -2.75; },
                                      [](I const& i, OBJ const& obj) { return false; }
    );
    auto constrain_p275 = ScalarObjectiveRankingParameter<E>(y.name(), OptimisationCriterion::MINIMISE, RankingConstraintSeverity::PERMISSIVE, y_p275,
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return ((o.evolve.bounding_box()[obj.variable].radius() - i.current_set.bounding_box()[obj.variable].radius())/(o.time-i.current_time)).get_d(); },
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return ((obj.radius - i.current_set.bounding_box()[obj.variable].radius())/(obj.time-i.current_time)).get_d(); },
                                      [](I const& i, OBJ const& obj) { return i.current_time > obj.time; }
    );
    auto constrain_m275 = ScalarObjectiveRankingParameter<E>(y.name(), OptimisationCriterion::MINIMISE, RankingConstraintSeverity::PERMISSIVE, y_m275,
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return ((o.evolve.bounding_box()[obj.variable].radius() - i.current_set.bounding_box()[obj.variable].radius())/(o.time-i.current_time)).get_d(); },
                                      [](I const& i, O const& o, DurationType const& d, OBJ const& obj) { return ((obj.radius - i.current_set.bounding_box()[obj.variable].radius())/(obj.time-i.current_time)).get_d(); },
                                      [](I const& i, OBJ const& obj) { return i.current_time > obj.time; }
    );

    auto step_size_used = ScalarRankingParameter<E>("step_size_used", OptimisationCriterion::MAXIMISE, [](I const& i,O const& o,DurationType const& d) { return o.step_size_used.get_d(); });

    List<Pair<TaskRankingParameter<E>,double>> specification = {{verification_p275,2}
                                                                ,{verification_m275,2}
                                                                ,{constrain_p275,2}
                                                                ,{constrain_m275,2}
                                                                /*,{execution_time_ranking<E>,1}
                                                                ,{step_size_used,2}*/
                                                                 };

    VerificationManager::instance().set_ranking_space(evolver,specification);

    Real x0 = 1.4_dec;
    Real y0 = 2.4_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    auto function_factory = TaylorFunctionFactory(sweeper1);
    EnclosureConfiguration enclosure_config(function_factory);
    enclosure_config.set_reconditioning_num_blocks(3);
    auto initial_set = evolver.enclosure({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0},enclosure_config);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,initial_set);

    auto unsafe_upper_set = evolver.enclosure({-2.5_dec<=x<=2.5_dec,-3<=y<=-2.75_dec},enclosure_config);
    auto unsafe_lower_set = evolver.enclosure({-2.5_dec<=x<=2.5_dec,2.75_dec<=y<=3},enclosure_config);

    Real evolution_time = 7;

    SizeType num_tries = 100;
    List<SizeType> success_times;
    List<SizeType> success_num_reachables;
    ProgressIndicator indicator(num_tries);
    List<List<TaskExecutionRanking>> list_best_rankings;
    auto start_all = std::chrono::high_resolution_clock::now();
    ARIADNE_LOG_SCOPE_CREATE;
    for (SizeType i=0; i<num_tries; ++i) {
        ARIADNE_LOG_PRINTLN_AT(1,"Test n." << i);
        auto start = std::chrono::high_resolution_clock::now();
        ARIADNE_LOG_PRINTLN_AT(1,"Computing orbit... ");
        try {
            auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
            auto end = std::chrono::high_resolution_clock::now();
            SizeType time_millis = (SizeType)std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            success_times.push_back(time_millis);
            success_num_reachables.push_back(orbit.reach().size());
            ARIADNE_LOG_PRINTLN_AT(1,"Done in " << ((double)time_millis)/1000 << " seconds.");

            ARIADNE_LOG_PRINTLN_AT(1,"Optimal point: " << ConcurrencyManager::instance().optimal_point());

            /*
            ARIADNE_LOG_PRINTLN("Plotting...");
            LabelledFigure fig({-2.5<=x<=2.5,-3<=y<=3});
            fig << fill_colour(1.0,0.75,0.5);
            fig.draw(orbit.reach());
            fig.draw(orbit.initial());
            fig.draw(unsafe_lower_set);
            fig.draw(unsafe_upper_set);
            fig.write("vanderpol");
             */

        } catch (CriticalRankingFailureException<VectorFieldEvolver>& ex) {
            ARIADNE_LOG_PRINTLN_AT(1,"Safety verification failure: " << ex.what());
        }
        list_best_rankings.push_back(ConcurrencyManager::instance().best_rankings());
        ConcurrencyManager::instance().clear_best_rankings();
        ConcurrencyManager::instance().choose_runner_for(evolver);
        VerificationManager::instance().set_ranking_space(evolver,specification);
        indicator.update_current(i+1);
        ARIADNE_LOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% of tries, currently " << 1+i-success_times.size() << " failures.");
    }
    double avg_time = 0; double stddev_time = 0;
    double avg_num_reachable = 0; double stddev_num_reachable = 0;
    for (SizeType i=0; i<success_times.size(); ++i) {
        auto time = success_times.at(i);
        auto num_reachable = success_num_reachables.at(i);
        avg_time += time; stddev_time += time*time;
        avg_num_reachable += num_reachable; stddev_num_reachable += num_reachable*num_reachable;
    }
    avg_time = avg_time/success_times.size()/1e3;
    stddev_time = sqrt(stddev_time/1e6/success_times.size() - avg_time*avg_time);
    avg_num_reachable = avg_num_reachable/success_num_reachables.size();
    stddev_num_reachable = sqrt(stddev_num_reachable/success_num_reachables.size() - avg_num_reachable*avg_num_reachable);

    ARIADNE_LOG_PRINTLN("Got " << round(((double)num_tries-success_times.size())/num_tries*1000)/10 <<
        "% of errors, with average execution time of " << round(avg_time*100)/100 << " s (±" << round(stddev_time*100)/100 << ")" <<
        " and average number of reachable sets " << round(avg_num_reachable) << " (±" << round(stddev_num_reachable) << ")");

    auto end_all = std::chrono::high_resolution_clock::now();
    SizeType time_millis = (SizeType)std::chrono::duration_cast<std::chrono::milliseconds>(end_all-start_all).count();
    ARIADNE_LOG_PRINTLN("Took a total of " << ((double)time_millis)/1000 << " seconds.");

    print_all_rankings(list_best_rankings);
}
