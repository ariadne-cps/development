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

    //RealConstant a("a",0.3_dec);
    RealVariable x("x"),y("y"),z("z");
    RealVariable a("a");
    RealVariable t("t");

    VectorField dynamics({dot(x)=-x*y/(x+1),
                          dot(y)=x*y/(x+1)-a*y,
                          dot(z)=a*y,
                          dot(a)=0,
                          dot(t)=1
                         });

    std::cout << "Production-destruction system:\n" << std::flush;

    MaximumError max_err=1e-2;
    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),max_err/1024/10);
    TaylorSeriesIntegrator integrator(max_err,sweeper,LipschitzConstant(0.5),Order(2u));

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.01);
    evolver.configuration().set_maximum_spacial_error(1e-2);
    evolver.verbosity = evolver_verbosity;

    Real eps = 0.00_dec;

    Box<RealInterval> initial_set({{10.0_dec,10.0_dec},{0.01_dec,0.01_dec},{0.01_dec,0.01_dec},{0.3_dec,0.3_dec},{0,0}});

    Real evolution_time(100.0);

    StopWatch sw;

    std::cout << "Computing orbit... " << std::endl << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);

    sw.click();
    std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    std::cout << "Plotting..." << std::endl;
    LabelledFigure fig_tx(Axes2d({0<=t<=100,-1<=x<=11}));
    fig_tx.draw(orbit.reach());
    fig_tx.write("production_destruction-tx");
    LabelledFigure fig_ty(Axes2d({0<=t<=100,-1<=y<=11}));
    fig_ty.draw(orbit.reach());
    fig_ty.write("production_destruction-ty");
    LabelledFigure fig_tz(Axes2d({0<=t<=100,-1<=z<=11}));
    fig_tz.draw(orbit.reach());
    fig_tz.write("production_destruction-tz");
    std::cout << "Files written." << std::endl;
}
