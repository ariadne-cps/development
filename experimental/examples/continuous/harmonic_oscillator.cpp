/***************************************************************************
 *            harmonic_oscillator.cpp
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

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));
    Logger::instance().configuration().set_theme(TT_THEME_DARK);
    Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);
    ConcurrencyManager::instance().set_concurrency(4);
    Logger::instance().use_blocking_scheduler();

    ARIADNE_LOG_PRINTLN("Harmonic Oscillator");

    RealConstant f("f",1);
    RealVariable x("x"), y("y");

    VectorField system({dot(x)=2*pi*f*y, dot(y)=-2*pi*f*x});

    double max_err = 1e-8;
    auto sweeper1 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/10);
    auto sweeper2 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/100);
    auto sweeper3 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err);
    auto integrator_configuration = Configuration<TaylorPicardIntegrator>()
            .set_step_maximum_error(1e-7,1e-4)
            .set_maximum_temporal_order(9,15)
            .set_lipschitz_tolerance(0.5)
            .set_sweeper(sweeper2);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,integrator_configuration.search_space());
    TaylorPicardIntegrator integrator(integrator_configuration);

    typedef VectorFieldEvolver E; typedef TaskInput<E> I; typedef TaskOutput<E> O; typedef TaskObjective<E> OBJ;

    E evolver(system,Configuration<E>().set_integrator(integrator));
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration());
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration().search_space());

    OBJ y_65(y,PositiveFloatDPUpperBound(FloatDP(0.25_x,DoublePrecision())),Dyadic(2));
    auto verification_parameter = ScalarRankingParameter<E>(y.name(), OptimisationCriterion::MINIMISE, [y](I const& i, O const& o, DurationType const& d) {
        return ((o.evolve.bounding_box()[y].radius()-i.current_set.bounding_box()[y].radius())/o.step_size_used).get_d(); });
    auto verification_constraint = TaskRankingConstraint<E>(verification_parameter, 1.25, RankingConstraintSeverity::CRITICAL);
    auto refinement = ConfigurationPropertyRefinement<E>(ConfigurationPropertyPath("integrator").append("step_maximum_error"),
                                                         {y_65},ProportionalRefiner(-1e-2));
    VerificationManager::instance().add_safety_specification(evolver, {verification_constraint}, {refinement});

    Real x0 = 0;
    Real y0 = 1;
    Real eps_x0 = 0.1_dec;
    Real eps_y0 = 0.1_dec;

    auto function_factory = TaylorFunctionFactory(sweeper1);
    EnclosureConfiguration enclosure_config(function_factory);
    enclosure_config.set_reconditioning_num_blocks(3);
    auto initial_set = evolver.enclosure({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0},enclosure_config);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,initial_set);

    Real evolution_time = 2.5_dec;

    auto start = std::chrono::high_resolution_clock::now();
    ARIADNE_LOG_PRINTLN("Computing orbit... ");
    try {
        auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
        auto end = std::chrono::high_resolution_clock::now();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << ((double)std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count())/1000 << " seconds.");

        ARIADNE_LOG_PRINTLN_AT(1,"Optimal point: " << ConcurrencyManager::instance().optimal_point());

        ARIADNE_LOG_PRINTLN("Plotting...");
        LabelledFigure fig({-2.0<=x<=2.0,-2<=y<=2});
        fig << fill_colour(1.0,0.75,0.5);
        fig.draw(orbit.reach());
        fig.write("harmonic-xy");
        LabelledFigure fig2({0<=TimeVariable()<=2.5_dec,-2.0<=y<=2.0});
        fig2 << fill_colour(1.0,0.75,0.5);
        fig2.draw(orbit.reach());
        fig2.write("harmonic-ty");
        LabelledFigure fig3({0<=TimeVariable()<=2.5_dec,-2.0<=x<=2.0});
        fig3 << fill_colour(1.0,0.75,0.5);
        fig3.draw(orbit.reach());
        fig3.write("harmonic-tx");
    } catch (CriticalRankingFailureException<E>& ex) {
        ARIADNE_LOG_PRINTLN("Safety verification failure: " << ex.what());
    }

    ConcurrencyManager::instance().print_best_rankings();
    ConcurrencyManager::instance().print_refinement_values();
}
