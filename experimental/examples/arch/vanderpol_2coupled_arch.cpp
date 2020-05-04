/***************************************************************************
 *            vanderpol_2coupled_arch.cpp
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

#include <cstdarg>
#include "ariadne.hpp"
#include "utility/stopwatch.hpp"

using namespace Ariadne;

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity=get_verbosity(argc,argv);

    RealVariable x1("x1"), y1("y1"), x2("x2"), y2("y2");

    std::cout << "Coupled van der Pol Oscillator system:\n" << std::flush;

    ListSet<LabelledEnclosure> reach1, reach2;
/*
    {
        std::cout << "Running for mu=1...\n" << std::flush;

        RealConstant mu("mu",1.0_dec);
        VectorField dynamics({dot(x1)=y1, dot(y1)=mu*(1-sqr(x1))*y1+x2-2*x1, dot(x2)=y2, dot(y2)=mu*(1-sqr(x2))*y2+x1-2*x2});

        MaximumError max_err = 1e-5;
        //TaylorSeriesIntegrator integrator(max_err, Order(5u));
        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(0.005);
        evolver.configuration().set_maximum_spacial_error(2e-4);
        evolver.verbosity = evolver_verbosity;

        Box<RealInterval> initial_set({{1.25_dec, 1.55_dec},{2.25_dec, 2.35_dec},{1.25_dec, 1.55_dec},{2.35_dec, 2.45_dec}});

        Real evolution_time(7.0);

        StopWatch sw;

        std::cout << "Computing orbit... \n" << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;

        std::cout << "Checking properties... \n" << std::flush;

        SizeType ce=0;
        for (auto set : orbit.reach()) {
            if (possibly(set.bounding_box().continuous_set()[1] >= 2.75_dec)) {
                std::cout << "set with y1=" << set.bounding_box().continuous_set()[1] << " is outside the specification." << std::endl;
                ++ce;
            }
            if (possibly(set.bounding_box().continuous_set()[3] >= 2.75_dec)) {
                std::cout << "set with y2=" << set.bounding_box().continuous_set()[3] << " is outside the specification." << std::endl;
                ++ce;
            }
        }
        sw.click();
        if (ce>0) std::cout << "Number of failures in satisfying the specification: " << ce << std::endl;
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;

        reach1.adjoin(orbit.reach());
    }
    */

    {
        std::cout << "Running for mu=2...\n" << std::flush;

        RealConstant mu("mu",2.0_dec);
        VectorField dynamics({dot(x1)=y1, dot(y1)=mu*(1-sqr(x1))*y1+x2-2*x1, dot(x2)=y2, dot(y2)=mu*(1-sqr(x2))*y2+x1-2*x2});

        MaximumError max_err = 1e-5;
        //TaylorSeriesIntegrator integrator(max_err, Order(4u));
        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.05);
        evolver.configuration().set_maximum_step_size(0.005);
        evolver.configuration().set_maximum_spacial_error(2e-4);
        evolver.verbosity = evolver_verbosity;

        Box<RealInterval> initial_set({{1.55_dec, 1.85_dec},{2.25_dec, 2.35_dec},{1.55_dec, 1.85_dec},{2.35_dec, 2.45_dec}});

        Real evolution_time(8.0);

        StopWatch sw;

        std::cout << "Computing orbit... \n" << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;

        std::cout << "Checking properties... \n" << std::flush;

        SizeType ce=0;
        for (auto set : orbit.reach()) {
            if (possibly(set.bounding_box().continuous_set()[1] >= 4.05_dec)) {
                std::cout << "set with y1=" << set.bounding_box().continuous_set()[1] << " is outside the specification." << std::endl;
                ++ce;
            }
            if (possibly(set.bounding_box().continuous_set()[3] >= 4.05_dec)) {
                std::cout << "set with y2=" << set.bounding_box().continuous_set()[3] << " is outside the specification." << std::endl;
                ++ce;
            }
        }
        sw.click();
        if (ce>0) std::cout << "Number of failures in satisfying the specification: " << ce << std::endl;
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;

        reach2.adjoin(orbit.reach());
    }


    std::cout << "Plotting..." << std::endl;

    Axes2d x1y1_axes(-2.5<=x1<=2.5,-4.05<=y1<=4.05);
    LabelledFigure fig1=LabelledFigure(x1y1_axes);
    fig1 << fill_colour(Colour(0.6,0.6,0.6));
    fig1.draw(reach1);
    fig1 << fill_colour(ariadneorange);
    fig1.draw(reach2);
    fig1.write("coupled-vanderpol-x1y1");
    Axes2d x2y2_axes(-2.5<=x2<=2.5,-4.05<=y2<=4.05);
    LabelledFigure fig2=LabelledFigure(x1y1_axes);
    fig2 << fill_colour(Colour(0.6,0.6,0.6));
    fig2.draw(reach1);
    fig2 << fill_colour(ariadneorange);
    fig2.draw(reach2);
    fig2.write("coupled-vanderpol-x2y2");
    std::cout << "Png files written." << std::endl;
}
