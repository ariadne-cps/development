/***************************************************************************
 *            henon_map.cpp
 *
 *  Copyright  2006-18  Pieter Collins
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
#include "utility/attribute.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression.hpp"
#include "function/taylor_model.hpp"
#include "algebra/differential.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "solvers/solver.hpp"
#include "geometry/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "dynamics/map.hpp"
#include "dynamics/map_evolver.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "geometry/grid_paving.hpp"
#include "dynamics/reachability_analyser.hpp"

    #include "geometry/function_set.hpp"
    #include "geometry/affine_set.hpp"

using namespace Ariadne;
using namespace std;

#define ARIADNE_PRINT(expr) { std::cout << #expr << ": " << (expr) << "\n"; }


Int main()
{

    // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
    Real a=Decimal(1.3), b=Decimal(0.3);
    RealVariable x("x"), y("y");
    EffectiveVectorMultivariateFunction henon = make_function({x,y},{a-x*x+b*y,x});
    ARIADNE_PRINT(henon);

    // Compute a fixed point
    IntervalNewtonSolver solver(maximum_error=1e-2, maximum_number_of_steps=16);
    ExactBoxType search_box({{0,1},{0,1}});
    Point<FloatDPBounds> fixed_point = Point(solver.fixed_point(henon,search_box));
    ARIADNE_PRINT(fixed_point);


    // Set up the evaluators
    MapEvolver evolver(henon);
    ReachabilityAnalyser<IteratedMap> analyser(henon,evolver);
    analyser.configuration().set_maximum_grid_depth(5);
    analyser.set_verbosity(3);

    // Set-up initial set and time for evolution
    RealBoxSet initial_set={{0.5_dec,0.6_dec},{0.95_dec,1.05_dec}};
    Enclosure initial_enclosure(initial_set,evolver.enclosure_configuration());

    // Set up the evolution parameters and grid
    IteratedMap::TimeType time(4);

    // Compute the reachable sets
    ListSet<Enclosure> evolve_set,reach_set;
    tie(reach_set,evolve_set) = evolver.reach_evolve(initial_enclosure,time);
    reach_set = evolver.reach(initial_enclosure,time);
    evolve_set = evolver.evolve(initial_enclosure,time);

    ExactBoxType figure_bounding_box({{-4,2},{-3,3}});
    Figure fig(figure_bounding_box,0,1);
    fig << fill_colour(magenta) << reach_set;
    fig << fill_colour(blue) << initial_set;
    fig << fill_colour(red) << evolve_set;
    fig << line_width(4) << fill_colour(green) << fixed_point << line_width(1);
    fig.write("henon_map-reach");

    // Compute the chain-reach set
    GridTreePaving chain_reach_set = analyser.outer_chain_reach(initial_set);

    fig.clear();
    fig << fill_colour(blue) << initial_set << fill_colour(cyan) << chain_reach_set;
    fig.write("henon_map-chain_reach");
}
