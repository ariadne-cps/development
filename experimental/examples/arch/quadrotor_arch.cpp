/***************************************************************************
 *            quadrotor.cpp
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

    RealConstant g("g",9.81_dec);
    RealConstant R("R",0.1_dec);
    RealConstant l("l",0.5_dec);
    RealConstant Mrotor("Mrotor",0.1_dec);
    RealConstant M("M",1);
    RealConstant m("m",M+4*Mrotor);
    RealConstant Jx("Jx",2*M*sqr(R)/5 + 2*sqr(l)*Mrotor);
    RealConstant Jy("Jy",Jx);
    RealConstant Jz("Jz",2*M*sqr(R)/5 + 4*sqr(l)*Mrotor);

    RealVariable x1("x1"),x2("x2"),x3("x3"),x4("x4"),x5("x5"),x6("x6"),x7("x7"),x8("x8"),x9("x9"),x10("x10"),x11("x11"),t("t");

    VectorField dynamics({dot(x1)=cos(x8)*cos(x9)*x4+(sin(x7)*sin(x8)*cos(x9)-cos(x7)*sin(x9))*x5+(cos(x7)*sin(x8)*cos(x9)+sin(x7)*sin(x9))*x6,
                          dot(x2)=cos(x8)*sin(x9)*x4+(sin(x7)*sin(x8)*sin(x9)+cos(x7)*cos(x9))*x5+(cos(x7)*sin(x8)*sin(x9)-sin(x7)*cos(x9))*x6,
                          dot(x3)=sin(x8)*x4-sin(x7)*cos(x8)*x5-cos(x7)*cos(x8)*x6,
                          dot(x4)=x1*x5-x11*x6-g*sin(x8),
                          dot(x5)=x10*x6+g*cos(x8)*sin(x7),
                          dot(x6)=x11*x4-x10*x5+g*cos(x8)*cos(x7)-(m*g-10*(x3-1)+3*x6)/m,
                          dot(x7)=x10+sin(x7)*tan(x8)*x11,
                          dot(x8)=cos(x7)*x11,
                          dot(x9)=sin(x7)/cos(x8)*x11,
                          dot(x10)=1/Jx*(-x7-x10),
                          dot(x11)=1/Jy*(-x8-x11),
                          dot(t)=Real(1.0)
                         });

    std::cout << "Quadrotor system:\n" << std::flush;

    MaximumError max_err=1e-2;
    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.01);
    evolver.configuration().set_maximum_spacial_error(1e-2);
    evolver.verbosity = evolver_verbosity;

    Real evolution_time(5.0);

    ListSet<LabelledEnclosure> reach1, reach2, reach3;

    {
        Real eps = 0.1_dec;

        RealVariablesBox initial_set(
                {-eps <= x1 <= eps, -eps <= x2 <= eps, -eps <= x3 <= eps, -eps <= x4 <= eps, -eps <= x5 <= eps,
                 -eps <= x6 <= eps,
                 x6 == 0, x7 == 0, x8 == 0, x9 == 0, x10 == 0, x11 == 0, t == 0});

        StopWatch sw;

        std::cout << "Computing orbit for Delta=0.1 ... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach1 = orbit.reach();
        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    {
        Real eps = 0.4_dec;

        RealVariablesBox initial_set(
                {-eps <= x1 <= eps, -eps <= x2 <= eps, -eps <= x3 <= eps, -eps <= x4 <= eps, -eps <= x5 <= eps,
                 -eps <= x6 <= eps,
                 x6 == 0, x7 == 0, x8 == 0, x9 == 0, x10 == 0, x11 == 0, t == 0});

        StopWatch sw;

        std::cout << "Computing orbit for Delta=0.4 ... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach2 = orbit.reach();
        std::cout << "Checking properties... " << std::endl << std::flush;

        for (auto set : orbit.reach()) {
            if (possibly(set.bounding_box().continuous_set()[2] >= 1.40_dec))
                std::cout << "height of " << set.bounding_box().continuous_set()[2] << " is over the required bound."
                          << std::endl;
            if (possibly(set.bounding_box().continuous_set()[11] >= 1) and
                possibly(set.bounding_box().continuous_set()[2] <= 0.9_dec))
                std::cout << "height of " << set.bounding_box().continuous_set()[2]
                          << " is below the required bound after 1s." << std::endl;
            if (possibly(set.bounding_box().continuous_set()[11] >= 5) and possibly(
                    set.bounding_box().continuous_set()[2] <= 0.98_dec or
                    set.bounding_box().continuous_set()[2] >= 1.02_dec))
                std::cout << "height of " << set.bounding_box().continuous_set()[2]
                          << " is outside the required bounds at 5s." << std::endl;
        }
        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    {
        Real eps = 0.8_dec;

        RealVariablesBox initial_set(
                {-eps <= x1 <= eps, -eps <= x2 <= eps, -eps <= x3 <= eps, -eps <= x4 <= eps, -eps <= x5 <= eps,
                 -eps <= x6 <= eps,
                 x6 == 0, x7 == 0, x8 == 0, x9 == 0, x10 == 0, x11 == 0, t == 0});

        StopWatch sw;

        std::cout << "Computing orbit for Delta=0.8 ... " << std::endl << std::flush;
        auto orbit = evolver.orbit(evolver.enclosure(initial_set), evolution_time, Semantics::UPPER);
        reach3 = orbit.reach();
        sw.click();
        std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;
    }

    std::cout << "Plotting..." << std::endl;
    LabelledFigure fig(Axes2d({0<=t<=5,-0.8<=x3<=1.5}));
    fig << line_colour(0.0,0.0,0.0);
    fig << line_style(false);
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(reach3);
    fig << fill_colour(0.6,0.6,0.6);
    fig.draw(reach2);
    fig << fill_colour(1.0,1.0,1.0);
    fig.draw(reach1);
    fig.write("quadrotor");
    std::cout << "File quadrotor.png written." << std::endl;
}
