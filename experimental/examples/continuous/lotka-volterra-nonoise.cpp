/***************************************************************************
 *            lotka-volterra-nonoise.cpp
 *
 *  Copyright  2017  Luca Geretti
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

using namespace Ariadne;


int main(int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));

    RealVariable u1("u1"), u2("u2");
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=u1*x*(1-y), dot(y)= u2*y*(x-1), dot(u1)=0, dot(u2)=0});

    MaximumError max_err=1e-5;
    TaylorPicardIntegrator integrator(max_err);
    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(1.0/50);
    evolver.configuration().set_maximum_spacial_error(max_err);

    RealVariablesBox initial_set({x==1.2_dec,y==1.1_dec,2.99_dec<=u1<=3.01_dec,0.99_dec<=u2<=1.01_dec});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time = 10;

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    auto bbx = orbit.final().bounding_box();
    double volume = (bbx[x].width()*bbx[y].width()).get_d();
    std::cout << "Volume cost: " << 1.0/std::pow(volume,1.0/bbx.dimension()) << std::endl;

    Axes2d axes(0.7<=x<=1.4,0.7<=y<=1.4);
    LabelledFigure fig=LabelledFigure(axes);
    fig << line_colour(0.0,0.0,0.0);
    fig << line_style(true);
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(orbit.reach());
    fig.write("lotka-volterra");
}
