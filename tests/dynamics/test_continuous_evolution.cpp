/***************************************************************************
 *            test_continuous_evolution.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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

#include <fstream>
#include <iostream>

#include "config.hpp"
#include "utility/tuple.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "solvers/integrator.hpp"
#include "symbolic/expression_set.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/orbit.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "function/user_function.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

/// This function diverges heavily
struct FailOne : VectorMultivariateFunctionData<2,2,1> {
    template<class R, class A, class P> static Void
    compute(R& r, const A& x, const P& p) {
          r[0] = 1;
          r[1] = -p[0] * x[1] + p[0];
    }
};

/// This function diverges heavily
struct FailTwo : VectorMultivariateFunctionData<3,3,1> {
    template<class R, class A, class P> static Void
    compute(R& r, const A& x, const P& p) {
          r[0] = 1;
          r[1] = x[1] * x[2] / p[0];
          r[2] = 0;
    }
};


class TestContinuousEvolution
{
  public:
    Void test() const;
    Void failure_test() const;
};

Int main()
{
    //std::cerr<<"SKIPPED "; return 1;
    ARIADNE_TEST_CALL(TestContinuousEvolution().test());
    //ARIADNE_TEST_CALL(TestContinuousEvolution().failure_test());
    return ARIADNE_TEST_FAILURES;
}



Void TestContinuousEvolution::test() const
{
    // cout << __PRETTY_FUNCTION__ << endl;

    typedef VectorField::EnclosureType EnclosureType;

    // Set up the evolution parameters and grid
    Real time(2.0_dec);
    ExactDouble step_size(0.5_x);
    ExactDouble enclosure_radius(0.25_x);

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-8_pr);

    // Set up the evaluators

    Configuration<TaylorPicardIntegrator> integrator_configuration;
    integrator_configuration.set_step_maximum_error(1e-2);
    TaylorPicardIntegrator picard_integrator(Configuration<TaylorPicardIntegrator>()
        .set_step_maximum_error(1e-6)
        .set_sweeper(sweeper)
        .set_maximum_temporal_order(8)
    );

    // Set up the evaluators
    GradedTaylorSeriesIntegrator series_integrator(maximum_error=1e-4_pr,sweeper,lipschitz_constant=0.5_x,step_maximum_error=1e-6_pr,
                                             minimum_spacial_order=1,minimum_temporal_order=4,maximum_spacial_order=3,maximum_temporal_order=8);

    IntegratorInterface& integrator=picard_integrator;

    ARIADNE_TEST_PRINT(integrator);

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Real mu=Dyadic(0.5_x);
    RealVariable x("x"), v("v");

    VectorField vanderpol({dot(x)=v,dot(v)=mu*(1-x*x)*v-x});
    ARIADNE_TEST_PRINT(vanderpol);

    // Define the initial box
    RealVariablesBox initial_box({x.in(1.01_dec,1.02_dec),v.in(0.51_dec,0.52_dec)});
//    initial_box[0]=ExactIntervalType(1.01,1.02);
//    initial_box[1]=ExactIntervalType(0.51,0.52);

    VectorFieldEvolver evolver(vanderpol,Configuration<VectorFieldEvolver>()
        .set_integrator(integrator)
        .set_enable_reconditioning(true)
        .set_maximum_spacial_error(1e-2)
        .set_maximum_enclosure_radius(enclosure_radius)
        .set_maximum_step_size(step_size)
    );

    // Over-approximate the initial set by a grid cell
    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
    EnclosureType initial_set(initial_box,vanderpol.state_space(),EnclosureConfiguration(function_factory));
    ARIADNE_TEST_PRINT(initial_set);

    Semantics semantics=Semantics::UPPER;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolver.orbit(initial_set,time,semantics);
    ARIADNE_TEST_PRINT(orbit);

/*
    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    Figure fig;
    fig.set_bounding_box(ExactBoxType{{-1.0,+21.0},{-1.125,+1.125}});
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-vdp");
*/

//    LabelledFigure fig(Axes2d(-1.0<=x<=+21,-1.125<=v<=+1.125));
    LabelledFigure fig(Axes2d(-1.0<=x<=+21,-1.125<=v<=+1.125));
//    fig.set_bounds({x,{-1.0_dec,+21.0_dec}},{v,{-1.125_dec,+1.125_dec}});
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-vdp");
}

Void TestContinuousEvolution::failure_test() const
{
    // The systems in this test are stiff and are expected to fail.

    // cout << __PRETTY_FUNCTION__ << endl;

    typedef VectorField::EnclosureType EnclosureType;

    // Set up the evolution parameters and grid
    Real time(0.5_dec);
    ExactDouble step_size(0.01_pr);
    ExactDouble enclosure_radius(0.25_x);

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-8_pr);

    // Set up the evaluators
    TaylorPicardIntegrator integrator(Configuration<TaylorPicardIntegrator>()
        .set_step_maximum_error(1e-8)
        .set_sweeper(sweeper)
        .set_lipschitz_tolerance(0.5)
        .set_minimum_temporal_order(0)
        .set_maximum_temporal_order(6)
    );

    // Define the initial box
    ExactBoxType initial_box = ExactBoxType{{0.0_x,0.0_x},{0.9_pr,0.9_pr}};

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field for the first test
    Real p = 200;
    RealVariable x("x"), y("y");

    VectorField failone_vf({dot(x)=1,dot(y)=-p*y+p});
    VectorFieldEvolver evolverone(failone_vf,Configuration<VectorFieldEvolver>()
        .set_integrator(integrator)
        .set_enable_reconditioning(true)
        .set_maximum_spacial_error(1e-6)
        .set_maximum_enclosure_radius(enclosure_radius)
        .set_maximum_step_size(step_size)
    );

    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-10));
    EnclosureType initial_set(initial_box,failone_vf.state_space(),EnclosureConfiguration(function_factory));
    // cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=Semantics::UPPER;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolverone.orbit(initial_set,time,semantics);

    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
    cout << "final set radius="<< orbit.final()[0].radius() << endl;

    ARIADNE_TEST_COMPARE(orbit.final()[0].radius(),>,0.5);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..."<< std::flush;
    LabelledFigure fig(Variables2d(x,y),widen(orbit.reach().bounding_box(),+0.25_x));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failone");
/*
    Figure fig; fig.set_bounding_box( widen(orbit.reach().bounding_box(),+0.25_x));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach().euclidean_set();
    fig << fill_colour(magenta) << orbit.intermediate().euclidean_set();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failone");
*/
    // Set up the vector field for the second test
    p = Rational(1,10);
    EffectiveScalarMultivariateFunction z3=EffectiveScalarMultivariateFunction::zero(3);
    EffectiveScalarMultivariateFunction o3=EffectiveScalarMultivariateFunction::constant(3,1.0_q);
    EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
    EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
    EffectiveScalarMultivariateFunction x2=EffectiveScalarMultivariateFunction::coordinate(3,1);
    EffectiveVectorMultivariateFunction failtwo={o3,x1*x2/p,z3};
    VectorField failtwo_vf(failtwo);
    VectorFieldEvolver evolvertwo(failtwo_vf,Configuration<VectorFieldEvolver>()
        .set_integrator(integrator)
        .set_enable_reconditioning(true)
        .set_maximum_spacial_error(1e-2)
        .set_maximum_enclosure_radius(enclosure_radius)
        .set_maximum_step_size(step_size)
    );

    ExactBoxType initial_box2 = ExactBoxType{{0.0_x,0.0_x},{1.0_x,1.0_x},{1.0_x,1.0_x}};
    initial_set = EnclosureType(initial_box2,failtwo_vf.state_space(),EnclosureConfiguration(function_factory));

    time = 1.5_dec;

    // Compute the reachable sets
    orbit = evolvertwo.orbit(initial_set,time,semantics);

    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
    cout << "final set radius="<< orbit.final()[0].radius() << endl;

    ARIADNE_TEST_COMPARE(orbit.final()[0].radius(),>,1.0);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    fig.clear(); fig.set_bounding_box(orbit.reach().bounding_box());
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failtwo");
    std::cout << std::endl;
}

