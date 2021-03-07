/***************************************************************************
 *            higgins-selkov-nonoise.cpp
 *
 *  Copyright  2018  Luca Geretti
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

    RealVariable v0("v0"), k1("k1"), k2("k2");
    RealVariable S("S"), P("P");

    VectorField dynamics({dot(S)=v0-S*k1*P*P, dot(P)= S*k1*P*P-k2*P, dot(v0)=0, dot(k1)=0, dot(k2)=0});

    MaximumError max_err=1e-5;
    TaylorPicardIntegrator integrator(max_err);
    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(1.0/50);
    evolver.configuration().set_maximum_spacial_error(max_err);

    Real S0 = 2, P0 = 1, v00 = 1, k10 = 1, k20 = 1.00001_dec;
    Real eS = 0.01_dec, eP = 0.01_dec, en = 0.0002_dec;
    RealVariablesBox initial_set({S0-eS<=S<=S0+eS,P0-eP<=P<=P0+eP,v00-en<=v0<=v00+en,k10-en<=k1<=k10+en,k20-en<=k2<=k20+en});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time = 10;
    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    auto bbx = orbit.final().bounding_box();
    double volume = (bbx[S].width()*bbx[P].width()).get_d();
    std::cout << "Volume cost: " << 1.0/std::pow(volume,1.0/bbx.dimension()) << std::endl;

    Axes2d axes(0.5<=S<=1.5,0.5<=P<=1.5);
    LabelledFigure fig=LabelledFigure(axes);
    fig << line_colour(0.0,0.0,0.0);
    fig << line_style(true);
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(orbit.reach());
    fig.write("higgins-selkov");
}
