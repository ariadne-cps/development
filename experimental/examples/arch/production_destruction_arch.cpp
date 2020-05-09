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

Int main(Int argc, const char* argv[]) {
    Nat evolver_verbosity = get_verbosity(argc, argv);

    RealVariable x("x"), y("y"), z("z");
    RealVariable a("a");
    RealVariable t("t");

    VectorField dynamics({dot(x) = -x * y / (x + 1),
                                 dot(y) = x * y / (x + 1) - a * y,
                                 dot(z) = a * y,
                                 dot(a) = 0,
                                 dot(t) = 1
                         });

    std::cout << "Production-destruction system:\n" << std::flush;

    MaximumError max_err = 1e-2;
    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(), max_err / 1024 / 10);
    TaylorSeriesIntegrator integrator(max_err, sweeper, LipschitzConstant(0.5), Order(2u));
    //TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics, integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.01);
    evolver.configuration().set_maximum_spacial_error(1e-2);
    evolver.verbosity = evolver_verbosity;

    ListSet<Enclosure> reach1, reach2, reach3;

    Real evolution_time(100.0);

    {
        Box<RealInterval> initial_set({{9.5_dec, 10.0_dec},
                                       {0.01_dec, 0.01_dec},
                                       {0.01_dec, 0.01_dec},
                                       {0.3_dec,  0.3_dec},
                                       {0,        0}});
        StopWatch sw;

        std::cout << "Computing orbit for 'I' setup... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach1 = orbit.reach();
        std::cout << "Done." << std::endl;

        Enclosure final = orbit.final()[0];
        auto bbox = final.bounding_box();
        if (bbox[0].midpoint().get_d() < 0.0)
            std::cout << "x is not >= 0" << std::endl;
        if (bbox[1].midpoint().get_d() < 0.0)
            std::cout << "y is not >= 0" << std::endl;
        if (bbox[2].midpoint().get_d() < 0.0)
            std::cout << "z is not >= 0" << std::endl;

        auto sum = bbox[0]+bbox[1]+bbox[2];
        if (definitely(not contains(sum,100.0)))
            std::cout << "x+y+z does not contain 100" << std::endl;

        auto volume = bbox[0].width()*bbox[1].width()*bbox[2].width();
        std::cout << "volume = " << volume << std::endl;

        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    {
        Box<RealInterval> initial_set({{10.0_dec, 10.0_dec},
                                       {0.01_dec, 0.01_dec},
                                       {0.01_dec, 0.01_dec},
                                       {0.296_dec, 0.304_dec},
                                       {0,        0}});
        StopWatch sw;

        std::cout << "Computing orbit for 'P' setup... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach2 = orbit.reach();
        std::cout << "Done." << std::endl;

        Enclosure final = orbit.final()[0];
        auto bbox = final.bounding_box();
        if (bbox[0].midpoint().get_d() < 0.0)
            std::cout << "x is not >= 0" << std::endl;
        if (bbox[1].midpoint().get_d() < 0.0)
            std::cout << "y is not >= 0" << std::endl;
        if (bbox[2].midpoint().get_d() < 0.0)
            std::cout << "z is not >= 0" << std::endl;

        auto sum = bbox[0]+bbox[1]+bbox[2];
        if (definitely(not contains(sum,100.0)))
            std::cout << "x+y+z does not contain 100" << std::endl;

        auto volume = bbox[0].width()*bbox[1].width()*bbox[2].width();
        std::cout << "volume = " << volume << std::endl;

        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    {
        Box<RealInterval> initial_set({{9.7_dec, 10.0_dec},
                                       {0.01_dec, 0.01_dec},
                                       {0.01_dec, 0.01_dec},
                                       {0.298_dec, 0.302_dec},
                                       {0,        0}});
        StopWatch sw;

        std::cout << "Computing orbit for 'I+P' setup... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach3 = orbit.reach();
        std::cout << "Done." << std::endl;

        Enclosure final = orbit.final()[0];
        auto bbox = final.bounding_box();
        if (bbox[0].midpoint().get_d() < 0.0)
            std::cout << "x is not >= 0" << std::endl;
        if (bbox[1].midpoint().get_d() < 0.0)
            std::cout << "y is not >= 0" << std::endl;
        if (bbox[2].midpoint().get_d() < 0.0)
            std::cout << "z is not >= 0" << std::endl;

        auto sum = bbox[0]+bbox[1]+bbox[2];
        if (definitely(not contains(sum,100.0)))
            std::cout << "x+y+z does not contain 100" << std::endl;

        auto volume = bbox[0].width()*bbox[1].width()*bbox[2].width();
        std::cout << "volume = " << volume << std::endl;

        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    std::cout << "Plotting..." << std::endl;
    Box<FloatDPUpperInterval> graphics_box{{0.0,11.0},{0.0,11.0},{0.0,11.0},{0.0,0.3},{0.0,100.0}};
    Figure fig=Figure();
    fig.set_projection_map(Projection2d(5,4,2));
    fig.set_bounding_box(graphics_box);
    fig.set_line_colour(0.0,0.0,0.0);
    fig.set_line_style(false);
    fig.set_fill_colour(0.6,0.6,0.6);
    fig.draw(reach3);
    fig.set_fill_colour(1.0,1.0,1.0);
    fig.draw(reach1);
    fig.set_fill_colour(1.0,0.75,0.5);
    fig.draw(reach2);
    fig.write("production_destruction");
    std::cout << "File written." << std::endl;
}
