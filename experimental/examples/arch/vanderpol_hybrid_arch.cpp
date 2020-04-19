/***************************************************************************
 *            vanderpol_hybrid_arch.cpp
 *
 *  Copyright  2008-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <cstdarg>
#include "ariadne.hpp"
#include "utility/stopwatch.hpp"

using namespace Ariadne;
using std::cout; using std::endl; using std::flush;

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity=get_verbosity(argc,argv);

    /// Set the system parameters
    Real a = 0.5_dec;  // Coefficient of restitution
    Real g = 9.8_dec;

    /// Set the position and velocity functions.
    RealVariable x("x");
    RealVariable y("y");
    RealVariable cnt("cnt");

    RealConstant mu("mu",1.0_dec);
    RealConstant sqradius("sqradius",3.0_dec);

    StringVariable vanderpol("vanderpol");
    StringConstant outside("outside");
    StringConstant inside("inside");
    HybridAutomaton automaton(vanderpol.name());

    DiscreteEvent enter("enter");
    DiscreteEvent exit("exit");

    automaton.new_mode(vanderpol|outside,{dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x, dot(cnt)=0});
    automaton.new_mode(vanderpol|inside,{dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x, dot(cnt)=1});
    automaton.new_transition(vanderpol|outside,enter,vanderpol|inside,{next(x)=x,next(y)=y,next(cnt)=cnt},sqr(x)+sqr(y)<=sqradius,EventKind::IMPACT);
    automaton.new_transition(vanderpol|inside,exit,vanderpol|outside,{next(x)=x,next(y)=y,next(cnt)=cnt},sqr(x)+sqr(y)>=sqradius,EventKind::IMPACT);

    MaximumError max_err=1e-5;
    TaylorSeriesIntegrator integrator(max_err,Order(5u));

    GeneralHybridEvolver evolver(automaton);
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(2e-4);
    evolver.verbosity=evolver_verbosity;

    Real ex(0.01_dec); //Real ex(0.15_dec);
    Real ey(0.01_dec); //Real ey(0.05_dec);
    HybridSet initial_set(vanderpol|outside,{1.4_dec-ex<=x<=1.4_dec+ex,-ey+2.4_dec<=y<=2.4_dec+ey,cnt==0});
    HybridTime evolution_time(7.0,5);

    StopWatch sw;

    std::cout << "Computing evolution... " << std::flush;
    auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::LOWER);
    sw.click();
    std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;

    HybridAutomaton circle;
    DiscreteLocation rotate;
    RealVariable t("t");
    circle.new_mode(rotate,{let(x)=sqradius*cos(t),let(y)=-sqradius*sin(t)},{dot(t)=1});
    HybridSimulator simulator;
    simulator.set_step_size(0.1);
    HybridRealPoint circle_initial(rotate,{t=0});
    HybridTime circle_time(6.28_dec,1);
    auto circle_orbit = simulator.orbit(circle,circle_initial,circle_time);

    plot("vanderpol-xy",Axes2d(-2.5<=x<=2.5,-4<=y<=4), ariadneorange, orbit, black, circle_orbit);
    plot("vanderpol-tcnt",Axes2d(0<=TimeVariable()<=7,0<=cnt<=7), orbit);
}
