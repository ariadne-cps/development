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

    std::cout << "2 coupled van der Pol oscillators system:\n" << std::flush;

    ListSet<Enclosure> reach1, reach2;

    {
        std::cout << "Running for mu=1...\n" << std::flush;

        RealConstant mu("mu",1.0_dec);
        VectorField dynamics({dot(x1)=y1, dot(y1)=mu*(1-sqr(x1))*y1+x2-2*x1, dot(x2)=y2, dot(y2)=mu*(1-sqr(x2))*y2+x1-2*x2});

        MaximumError max_err = 1e-5;
        TaylorSeriesIntegrator integrator(max_err, Order(4u));
        //TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics, integrator);
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(0.01);
        evolver.configuration().set_maximum_spacial_error(2e-4);
        evolver.verbosity = evolver_verbosity;

        Box<RealInterval> initial_set({{1.25_dec, 1.55_dec},{2.25_dec, 2.35_dec},{1.25_dec, 1.55_dec},{2.25_dec, 2.35_dec}});

        Real evolution_time(7.0);

        StopWatch sw;

        std::cout << "Computing orbit... \n" << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
/*
        std::cout << "Checking properties... \n" << std::flush;

        SizeType ce=0;
        for (auto set : orbit.reach()) {
            if (possibly(set.bounding_box()[1] >= 2.75_dec)) {
                std::cout << "set with y=" << set.bounding_box()[1] << " is outside the specification." << std::endl;
                ++ce;
            }
        }
        sw.click();
        if (ce>0) std::cout << "Number of sets not satisfying the specification: " << ce << std::endl;
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
*/
        reach1.adjoin(orbit.reach());
    }

    /*
    std::cout << "Plotting..." << std::endl;
    Box<FloatDPUpperInterval> graphics_box(4);
    graphics_box[0] = FloatDPUpperInterval(-4.0,4.0);
    graphics_box[1] = FloatDPUpperInterval(-4.0,4.0);
    graphics_box[2] = FloatDPUpperInterval(-4.0,4.0);
    graphics_box[3] = FloatDPUpperInterval(-4.0,4.0);
    Figure fig=Figure();
    fig.set_bounding_box(graphics_box);
    fig.set_projection_map(PlanarProjectionMap(4,0,1))
    fig.set_line_colour(black);
    fig.set_line_style(true);
    fig.set_fill_colour(0.6,0.6,0.6);
    fig.draw(reach2);
    fig.set_fill_colour(ariadneorange);
    fig.draw(reach1);
    fig.write("vanderpol-2coupled");
    */
    plot("vanderpol-2coupled-x1y1",PlanarProjectionMap(4,0,1),ApproximateBoxType({{-4.0,4.0},{-4.0,4.0},{-4.0,4.0},{-4.0,4.0}}), Colour(ariadneorange), reach1);
    plot("vanderpol-2coupled-x2y2",PlanarProjectionMap(4,2,3),ApproximateBoxType({{-4.0,4.0},{-4.0,4.0},{-4.0,4.0},{-4.0,4.0}}), Colour(ariadneorange), reach1);
    std::cout << "Png files written." << std::endl;
}
