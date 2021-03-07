/***************************************************************************
 *            test_discrete_evolution.cpp
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
#include "function/taylor_model.hpp"
#include "algebra/differential.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "dynamics/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "dynamics/map.hpp"
#include "dynamics/map_evolver.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestMapEvolver
{
  public:
    Void test() const;
};

Int main()
{
    ARIADNE_TEST_CALL(TestMapEvolver().test());
    return ARIADNE_TEST_FAILURES;
}

Void TestMapEvolver::test() const
{
    typedef DoublePrecision PR;
    DoublePrecision pr;

    typedef IteratedMap::EnclosureType EnclosureType;

    // Define the initial box
    ExactBoxType initial_box(2);
    initial_box[0]=ExactIntervalType(1.01_pr,1.03_pr);
    initial_box[1]=ExactIntervalType(0.51_pr,0.53_pr);

    ARIADNE_TEST_PRINT(initial_box);

    // Set up the map field
    // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
    Real a=1.5_dyadic; Real b=0.375_dyadic;
    RealVariable x("x"), y("y");
    IteratedMap henon({ next(x)=a-x*x+b*y, next(y)=x });
    ARIADNE_TEST_PRINT(henon);

    // Function evaluation sanity check
    Vector<FloatApproximation<PR>> p={{a,b},pr};
    Vector<FloatApproximation<PR>> xa={{0.5_x,0.25_x},pr};
    Vector<FloatApproximation<PR>> hxa={p[0]-xa[0]*xa[0]+xa[1]*p[1], xa[0]};
    ARIADNE_TEST_EQUAL(henon.update_function().evaluate(xa),hxa);
    Matrix<FloatApproximation<PR>> dhxa={{-2*xa[0],p[1]},{{1.0,pr},{0.0,pr}}};
    ARIADNE_TEST_EQUAL(henon.update_function().jacobian(xa),dhxa);


    // Function evaluation sanity check
    ARIADNE_TEST_PRINT(initial_box);
    ARIADNE_TEST_PRINT(image(initial_box,henon.update_function()));
    ARIADNE_TEST_PRINT(jacobian_range(henon.update_function(),cast_vector(initial_box)));



    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box,henon.state_space(),EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-10))));
    ARIADNE_TEST_PRINT(initial_set);

    // Set up the evolution parameters and grid
    IteratedMap::TimeType time(3);
    double enclosure_radius(0.25);

    // Set up the evaluators
    MapEvolver evolver(henon);
    evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);

    // Compute the reachable sets
    ListSet<EnclosureType> evolve_set,reach_set;
    //evolve_set = evolver.evolve(initial_set,time);
    reach_set = evolver.reach(initial_set,time);
    ARIADNE_TEST_PRINT(initial_set.bounding_box());
    //cout << "evolve_bounding_boxes=" << evolve_set.bounding_boxes() << endl;
    ARIADNE_TEST_PRINT(reach_set.bounding_boxes());

/*
    // Print the intial, evolve and reach sets
    Figure fig;
    fig.set_bounding_box({{-4,2},{-3,3}});
    fig << line_style(true) << fill_colour(cyan) << reach_set;
    fig << fill_colour(yellow) << evolve_set;
    fig << fill_colour(blue) << initial_set;
    fig.write("test_discrete_evolution-henon");
*/

    // Print the intial, evolve and reach sets
    LabelledFigure lfig(Axes2d(-4,x,2, -3,y,3));
    lfig << line_style(true) << fill_colour(cyan) << reach_set;
    lfig << fill_colour(yellow) << evolve_set;
    lfig << fill_colour(blue) << initial_set;
    lfig.write("test_discrete_evolution-henon-labelled");
}
