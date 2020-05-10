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

    MaximumError max_err = 1e-6;
    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics, integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.04);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    evolver.verbosity = evolver_verbosity;

    ListSet<LabelledEnclosure> reach1, reach2, reach3;

    Real evolution_time(100.0);

    {
        RealVariablesBox initial_set({9.5_dec<=x<=10,y==0.01_dec,z==0.01_dec,a==0.3_dec,t==0});

        StopWatch sw;

        std::cout << "Computing orbit for 'I' setup... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach1 = orbit.reach();
        std::cout << "Done." << std::endl;

        auto final = orbit.final()[0];
        auto bbox = final.bounding_box();
        if (bbox.continuous_set()[0].midpoint().get_d() < 0.0)
            std::cout << "x is not >= 0" << std::endl;
        if (bbox.continuous_set()[1].midpoint().get_d() < 0.0)
            std::cout << "y is not >= 0" << std::endl;
        if (bbox.continuous_set()[2].midpoint().get_d() < 0.0)
            std::cout << "z is not >= 0" << std::endl;

        auto sum = bbox.continuous_set()[0]+bbox.continuous_set()[1]+bbox.continuous_set()[2];
        if (definitely(not contains(sum,100.0)))
            std::cout << "x+y+z does not contain 100" << std::endl;

        auto volume = bbox.continuous_set()[0].width()*bbox.continuous_set()[1].width()*bbox.continuous_set()[2].width();
        std::cout << "Final volume = " << volume << std::endl;
	std::cout << "Range for final z = " << bbox.continuous_set()[2] << std::endl;    
	
	sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    {
        RealVariablesBox initial_set({x==10,y==0.01_dec,z==0.01_dec,0.296_dec<=a<=0.304_dec,t==0});

        StopWatch sw;

        std::cout << "Computing orbit for 'P' setup... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach2 = orbit.reach();
        std::cout << "Done." << std::endl;

        auto final = orbit.final()[0];
        auto bbox = final.bounding_box();
        if (bbox.continuous_set()[0].midpoint().get_d() < 0.0)
            std::cout << "x is not >= 0" << std::endl;
        if (bbox.continuous_set()[1].midpoint().get_d() < 0.0)
            std::cout << "y is not >= 0" << std::endl;
        if (bbox.continuous_set()[2].midpoint().get_d() < 0.0)
            std::cout << "z is not >= 0" << std::endl;

        auto sum = bbox.continuous_set()[0]+bbox.continuous_set()[1]+bbox.continuous_set()[2];
        if (definitely(not contains(sum,100.0)))
            std::cout << "x+y+z does not contain 100" << std::endl;

        auto volume = bbox.continuous_set()[0].width()*bbox.continuous_set()[1].width()*bbox.continuous_set()[2].width();
        std::cout << "Final volume = " << volume << std::endl;
	std::cout << "Range for final z = " << bbox.continuous_set()[2] << std::endl;    

        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    {
        RealVariablesBox initial_set({9.7_dec<=x<=10,y==0.01_dec,z==0.01_dec,0.298_dec<=a<=0.302_dec,t==0});

        StopWatch sw;

        std::cout << "Computing orbit for 'I+P' setup... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach3 = orbit.reach();
        std::cout << "Done." << std::endl;

        auto final = orbit.final()[0];
        auto bbox = final.bounding_box();
        if (bbox.continuous_set()[0].midpoint().get_d() < 0.0)
            std::cout << "x is not >= 0" << std::endl;
        if (bbox.continuous_set()[1].midpoint().get_d() < 0.0)
            std::cout << "y is not >= 0" << std::endl;
        if (bbox.continuous_set()[2].midpoint().get_d() < 0.0)
            std::cout << "z is not >= 0" << std::endl;

        auto sum = bbox.continuous_set()[0]+bbox.continuous_set()[1]+bbox.continuous_set()[2];
        if (definitely(not contains(sum,100.0)))
            std::cout << "x+y+z does not contain 100" << std::endl;

        auto volume = bbox.continuous_set()[0].width()*bbox.continuous_set()[1].width()*bbox.continuous_set()[2].width();
        std::cout << "Final volume = " << volume << std::endl;
	std::cout << "Range for final z = " << bbox.continuous_set()[2] << std::endl;    

        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    std::cout << "Plotting..." << std::endl;
    LabelledFigure fig(Axes2d({0<=t<=100,0<=z<=11}));
    fig << line_colour(0.0,0.0,0.0);
    fig << line_style(false);
    fig << fill_colour(0.6,0.6,0.6);
    fig.draw(reach2);
    fig << fill_colour(1.0,0.0,0.0);
    fig.draw(reach3);
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(reach1);
    fig.write("production_destruction");
    std::cout << "File production_destruction.png written." << std::endl;
}
