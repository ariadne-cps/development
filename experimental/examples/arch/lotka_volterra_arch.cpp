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

    RealConstant alpha("alpha",3);
    RealConstant beta("beta",3);
    RealConstant gamma("gamma",1);
    RealConstant delta("delta",1);
    RealConstant radius("radius",0.2_dec);

    StringVariable lotkavolterra("lotkavolterra");
    StringConstant outside("outside");
    StringConstant inside("inside");
    HybridAutomaton automaton(lotkavolterra.name());

    DiscreteEvent enter("enter");
    DiscreteEvent exit("exit");

    Real cx = gamma/delta;
    Real cy = alpha/beta;

    automaton.new_mode(lotkavolterra|outside,{dot(x)=alpha*x-beta*x*y, dot(y)= delta*x*y-gamma*y, dot(cnt)=0});
    automaton.new_mode(lotkavolterra|inside,{dot(x)=alpha*x-beta*x*y, dot(y)= delta*x*y-gamma*y, dot(cnt)=1});
    //automaton.new_transition(lotkavolterra|outside,enter,lotkavolterra|inside,{next(x)=x,next(y)=y,next(cnt)=cnt},sqr(x-cx)+sqr(y-cy)<=sqr(radius),EventKind::IMPACT);
    //automaton.new_transition(lotkavolterra|inside,exit,lotkavolterra|outside,{next(x)=x,next(y)=y,next(cnt)=cnt},sqr(x-cx)+sqr(y-cy)>=sqr(radius),EventKind::IMPACT);

    MaximumError max_err=1e-5;
    TaylorSeriesIntegrator integrator(max_err,Order(5u));

    GeneralHybridEvolver evolver(automaton);
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(2e-4);
    evolver.verbosity=evolver_verbosity;

    RealPoint ic({1.2_dec,1.1_dec});
    Real ex(0.0_dec); //Real ex(0.15_dec);
    Real ey(0.0_dec); //Real ey(0.05_dec);
    HybridSet initial_set(lotkavolterra|outside,{ic[0]-ex<=x<=ic[0]+ex,-ey+ic[1]<=y<=ic[1]+ey,cnt==0});
    HybridTime evolution_time(3.64,5);

    StopWatch sw;

    std::cout << "Computing evolution... " << std::flush;
    auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::LOWER);
    sw.click();
    std::cout << "Done in " << sw.elapsed() << " seconds." << std::endl;

    HybridAutomaton circle;
    DiscreteLocation rotate;
    RealVariable t("t");
    circle.new_mode(rotate,{let(x)=cx+radius*cos(t),let(y)=cy-radius*sin(t)},{dot(t)=1});
    HybridSimulator simulator;
    simulator.set_step_size(0.1);
    HybridRealPoint circle_initial(rotate,{t=0});
    HybridTime circle_time(6.28_dec,1);
    auto circle_orbit = simulator.orbit(circle,circle_initial,circle_time);

    plot("lotkavolterra-xy",Axes2d(0.75<=x<=1.3,0.75<=y<=1.3), ariadneorange, orbit, black, circle_orbit);
    plot("lotkavolterra-tcnt",Axes2d(0<=TimeVariable()<=evolution_time.continuous_time(),0<=cnt<=evolution_time.continuous_time()), orbit);
}
