/***************************************************************************
 *            production_destruction_arch.cpp
 *
 *  Copyright  2019  Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"
#include "utility/stopwatch.hpp"

using namespace Ariadne;

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity=get_verbosity(argc,argv);

    RealConstant a("a",0.3_dec);
    RealVariable x("x"),y("y"),z("z");

    VectorField dynamics({dot(x)=-x*y/(x+1),
                          dot(y)=x*y/(x+1)-a*y,
                          dot(z)=a*y
                         });

    std::cout << "Production-destruction system:\n" << std::flush;

    MaximumError max_err=1e-2;
    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),max_err/1024/10);
    TaylorSeriesIntegrator integrator(max_err,sweeper,LipschitzConstant(0.5),Order(2u));

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.1);
    evolver.configuration().set_maximum_spacial_error(1e-2);
    evolver.verbosity = evolver_verbosity;

    Real eps = 0.002_dec;

    Box<RealInterval> initial_set({{0.98_dec-eps,0.98_dec+eps},{0.01_dec-eps,0.01_dec+eps},{0.01_dec-eps,0.01_dec+eps}});

    Real evolution_time(100.0);

    StopWatch sw;

    std::cout << "Computing orbit... " << std::endl << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    sw.click();
    std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    std::cout << "Plotting..." << std::endl;
    plot("production_destruction-xy",PlanarProjectionMap(3,0,1),ApproximateBoxType({{0.0,1.0},{0.0,1.0},{0.0,1.0}}), Colour(1.0,0.75,0.5), orbit);
    plot("production_destruction-xz",PlanarProjectionMap(3,0,2),ApproximateBoxType({{0.0,1.0},{0.0,1.0},{0.0,1.0}}), Colour(1.0,0.75,0.5), orbit);
    plot("production_destruction-yz",PlanarProjectionMap(3,1,2),ApproximateBoxType({{0.0,1.0},{0.0,1.0},{0.0,1.0}}), Colour(1.0,0.75,0.5), orbit);
    std::cout << "Files written." << std::endl;
}
