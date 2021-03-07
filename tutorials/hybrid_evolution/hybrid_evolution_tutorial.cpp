/***************************************************************************
 *            hybrid_evolution_tutorial.cpp
 *
 *  Copyright  2020  Luca Geretti
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

#include <ariadne.hpp>

using namespace Ariadne;

//! [get_tank]
HybridAutomaton get_tank()
{
    // Declare the system constants
    RealConstant alpha("alpha",0.02_decimal);
    RealConstant beta("beta",0.3_decimal);

    // Declare the variables for the dynamics
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Create the tank automaton
    HybridAutomaton automaton("tank");

    // Declare a trivial discrete location (we use an empty label since there is only one location).
    DiscreteLocation flow;

    // The water level is always given by the same dynamic.
    // The inflow is controlled by the valve aperture, the outflow depends on the
    // pressure, which is proportional to the water height.
    automaton.new_mode(flow,{dot(height)=beta*aperture-alpha*height});

    return automaton;
}
//! [get_tank]

//! [get_valve]
HybridAutomaton get_valve()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4);

    // Declare the shared system variables
    RealVariable aperture("aperture");

    // Declare the events we use
    DiscreteEvent stop_opening("stop_opening");
    DiscreteEvent stop_closing("stop_closing");
    DiscreteEvent can_open("can_open");
    DiscreteEvent can_close("can_close");

    // Declare the variable for the automaton name
    StringVariable valve("valve");

    // Create the valve automaton
    HybridAutomaton automaton(valve.name());

    // Declare the values the valve variable can have
    DiscreteLocation opening(valve|"opening");
    DiscreteLocation closed(valve|"closed");
    DiscreteLocation opened(valve|"opened");
    DiscreteLocation closing(valve|"closing");

    // Define the algebraic equations for the opened/closed locations.
    automaton.new_mode(opened,{let(aperture)=1});
    automaton.new_mode(closed,{let(aperture)=0});
    // Define the differential equations for the opening/closing locations.
    automaton.new_mode(opening,{dot(aperture)=+1/T});
    automaton.new_mode(closing,{dot(aperture)=-1/T});

    // Define the transitions: source location, event and target location;
    // then a mix of reset, guard and event kind can be present; if the event kind
    // is not specified, then also the guard can't be specified: this implicitly
    // means that the event is an input event for this automaton.
    automaton.new_transition(closed,can_open,opening,{next(aperture)=aperture});
    automaton.new_transition(opening,stop_opening,opened,aperture>=1,EventKind::URGENT);
    automaton.new_transition(opened,can_close,closing,{next(aperture)=aperture});
    automaton.new_transition(closing,stop_closing,closed,aperture<=0,EventKind::URGENT);

    return automaton;
}
//! [get_valve]

//! [get_controller]
HybridAutomaton get_controller()
{
    // Declare some constants
    RealConstant hmin("hmin",5.75_decimal);
    RealConstant hmax("hmax",7.75_decimal);
    RealConstant delta("delta",0.02_decimal);

    // Declare the shared system variables
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent can_open("can_open");
    DiscreteEvent can_close("can_close");
    DiscreteEvent must_open("must_open");
    DiscreteEvent must_close("must_close");

    // Declare the variable for the automaton name
    StringVariable controller("controller");

    // Create the controller automaton
    HybridAutomaton automaton(controller.name());

    // Declare the locations for the controller
    DiscreteLocation rising(controller|"rising");
    DiscreteLocation falling(controller|"falling");

    // Instantiate modes for each location with no dynamics
    automaton.new_mode(rising);
    automaton.new_mode(falling);

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    automaton.new_invariant(falling,height>=hmin-delta,must_open);
    automaton.new_invariant(rising,height<=hmax+delta,must_close);

    // Specify the transitions, starting from the source location, according to an event, to a target location;
    // Following those arguments you specify a guard and whether the event is permissive or urgent.
    automaton.new_transition(falling,can_open,rising,height<=hmin+delta,EventKind::PERMISSIVE);
    automaton.new_transition(rising,can_close,falling,height>=hmax-delta,EventKind::PERMISSIVE);

    return automaton;
}
//! [get_controller]

//! [get_system]
CompositeHybridAutomaton get_system()
{
    CompositeHybridAutomaton watertank_system("watertank",{get_tank(),get_valve(),get_controller()});
    return watertank_system;
}
//! [get_system]


//! [simulate_evolution]
Void simulate_evolution(const CompositeHybridAutomaton& system)
{
    ARIADNE_LOG_SCOPE_CREATE;
    // Re-introduce the shared system variables required for the initial set
    RealVariable aperture("aperture");
    RealVariable height("height");
    TimeVariable time;

    StringVariable valve("valve");
    StringVariable controller("controller");

    String opened("opened");
    String rising("rising");

    // Create a simulator object
    HybridSimulator simulator;
    simulator.set_step_size(0.01);

    // Set an initial point for the simulation
    HybridRealPoint initial_point({valve|opened,controller|rising},{height=7});

    // Define the termination: continuous time and maximum number of transitions
    HybridTerminationCriterion termination(30.0_x,5);

    // Compute a simulation trajectory
    ARIADNE_LOG_PRINTLN("Computing simulation trajectory...");
    auto orbit = simulator.orbit(system,initial_point,termination);

    // Plot the simulation trajectory using all different projections
    ARIADNE_LOG_PRINTLN("Plotting simulation trajectory...");
    plot("simulation_t-height",Axes2d(0<=time<=30,5<=height<=9),orbit);
    plot("simulation_t-aperture",Axes2d(0<=time<=30,-0.1<=aperture<=1.1),orbit);
    plot("simulation_height-aperture",Axes2d(5<=height<=9,-0.1<=aperture<=1.1),orbit);
    ARIADNE_LOG_PRINTLN("Done computing and plotting simulation trajectory!");
}
//! [simulate_evolution]

//! [create_evolver]
GeneralHybridEvolver create_evolver(const CompositeHybridAutomaton& system)
{
    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(system);

    // Set the evolver configuration
    evolver.configuration().set_maximum_enclosure_radius(3.0);
    evolver.configuration().set_maximum_step_size(0.25);

    ARIADNE_LOG_PRINTLN("Evolver configuration: " << evolver.configuration());

    return evolver;
}
//! End of [create_evolver]

//! [compute_evolution]
Void compute_evolution(const GeneralHybridEvolver& evolver)
{
    ARIADNE_LOG_SCOPE_CREATE;
    // Re-introduce the shared system variables required for the initial set
    RealVariable aperture("aperture");
    RealVariable height("height");
    TimeVariable time;

    StringVariable valve("valve");
    StringVariable controller("controller");

    String opened("opened");
    String rising("rising");

    // Define the initial set, by supplying the location as a list of locations for each composed automata, and
    // the continuous set as a list of variable assignments for each variable controlled on that location
    // (the assignment can be either a singleton value using the == symbol or an interval using the <= symbols)
    HybridBoundedConstraintSet initial_set({valve|opened,controller|rising},{6.9_decimal<=height<=7});
    // Define the termination: continuous time and maximum number of transitions
    HybridTerminationCriterion termination(30,5);
    // Compute the orbit using upper semantics
    ARIADNE_LOG_PRINTLN("Computing evolution flow tube...");
    auto orbit = evolver.orbit(initial_set,termination,Semantics::UPPER);

    // Plot the flow tube using two different projections
    ARIADNE_LOG_PRINTLN("Plotting evolution flow tube...");
    plot("finite_evolution_t-height",Axes2d(0<=time<=30,5<=height<=9),orbit);
    plot("finite_evolution_t-aperture",Axes2d(0<=time<=30,-0.1<=aperture<=1.1),orbit);
    plot("finite_evolution_height-aperture",Axes2d(5<=height<=9,-0.1<=aperture<=1.1),orbit);
    ARIADNE_LOG_PRINTLN("Done computing and plotting evolution flow tube!");
}
//! [compute_evolution]

//! [create_analyser]
HybridReachabilityAnalyser create_analyser(const GeneralHybridEvolver& evolver)
{
    ARIADNE_LOG_SCOPE_CREATE
    // Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);

    //  Set the analyser configuration
    analyser.configuration().set_maximum_grid_fineness(6);
    analyser.configuration().set_lock_to_grid_time(5);

    ARIADNE_LOG_PRINTLN("Analyser configuration: " << analyser.configuration());

    return analyser;
}
//! [create_analyser]

//! [compute_reachability]
Void compute_reachability(const HybridReachabilityAnalyser& analyser)
{
    ARIADNE_LOG_SCOPE_CREATE;
    // Re-introduce the shared system variables required for the initial set
    RealVariable aperture("aperture");
    RealVariable height("height");

    StringVariable valve("valve");
    StringVariable controller("controller");

    String opened("opened");
    String rising("rising");

    // Define the initial set and final time
    HybridBoundedConstraintSet initial_set({valve|opened,controller|rising},{6.9_decimal<=height<=7});
    ARIADNE_LOG_PRINTLN("initial_set: " << initial_set);
    HybridTime final_time(30.0_x,5);
    ARIADNE_LOG_PRINTLN("final_time: " << final_time);

    // Compute over-approximation to infinite-time reachable set using upper semantics.
    ARIADNE_LOG_PRINTLN("Computing upper reach set...");
    auto upper_reach = analyser.upper_reach(initial_set,final_time);
    ARIADNE_LOG_PRINTLN("Plotting upper reach set...");
    plot("upper_reach",Axes2d(5<=height<=9,-0.1<=aperture<=1.1),upper_reach);
    ARIADNE_LOG_PRINTLN("Done computing and plotting upper reach set!");

    // Compute over-approximation to infinite-time reachable set using upper semantics.
    ARIADNE_LOG_PRINTLN("Computing outer chain reach set...");
    auto outer_chain_reach = analyser.outer_chain_reach(initial_set);
    ARIADNE_LOG_PRINTLN("Plotting outer chain reach set...");
    plot("outer_chain_reach",Axes2d(5<=height<=9,-0.1<=aperture<=1.1),outer_chain_reach);
    ARIADNE_LOG_PRINTLN("Done computing and plotting outer chain reach set!");
}
//! [compute_reachability]

//! [main]
Int main(Int argc, const char* argv[])
{
    // Acquire the verbosity value from the command line
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));

    // Create the composed automaton
    CompositeHybridAutomaton watertank_system=get_system();

    // Print the system description on the command line
    ARIADNE_LOG_PRINTLN("System: " << watertank_system);

    // Compute an approximate simulation of the system evolution
    simulate_evolution(watertank_system);

    // Create an evolver object
    auto evolver = create_evolver(watertank_system);

    // Compute the system evolution
    compute_evolution(evolver);

    // Create an analyser object
    auto analyser = create_analyser(evolver);

    // Compute the system reachability
    compute_reachability(analyser);
}
//! [main]
