/***************************************************************************
 *            rectifier.cpp
 *
 *  Copyright  2009-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
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

using namespace Ariadne;
using std::cout; using std::endl; using std::flush;

/// Function for plotting the orbit and reachability set
template<class SET> Void plot(const char* filename, const Nat& xaxis, const Nat& yaxis, const Nat& numVariables, const ExactBoxType& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) {
    // Assigns local variables
    Figure fig;

    fig.set_projection(numVariables,xaxis,yaxis);
    fig.set_bounding_box(bbox);

    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
        // The rectangle to be drawn
        ApproximateBoxType rect (numVariables);
        // Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

        // Gets the number of times each variable interval would be divided by 2
        Nat numDivisions = static_cast<Nat>(MAX_GRID_DEPTH) / numVariables;
        // Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
        ApproximateDouble step_x = 1.0/(1 << (numDivisions + ((static_cast<Nat>(MAX_GRID_DEPTH) - numDivisions*numVariables > xaxis) ? 1 : 0)));
        // Initiates the x position to the bounding box left bound
        ApproximateDouble pos_x = bbox[0].lower_bound().get_d();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
        rect[yaxis] = bbox[1];
        // While between the interval
        while (pos_x.get_d() < bbox[0].upper_bound().get_d())
        {
            rect[xaxis] = ApproximateIntervalType(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
            pos_x += step_x; // Shifts the x position
            fig << rect; // Appends the rectangle
        }

        // Repeats for the rectangles in the y direction
        ApproximateDouble step_y = 1.0/(1 << (numDivisions + ((static_cast<Nat>(MAX_GRID_DEPTH) - numDivisions*numVariables > yaxis) ? 1 : 0)));
        ApproximateDouble pos_y = bbox[1].lower_bound().get_d();
        rect[xaxis] = bbox[0];
        while (pos_y.get_d() < bbox[1].upper_bound().get_d())
        {
            rect[yaxis] = ApproximateIntervalType(pos_y,pos_y+step_y);
            fig << rect;
            pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc);
    fig << set;
    fig.write(filename);
}

Int main(Int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));

    Real amplitude(4.0_dec);
    Real frequency(50.0_dec);
    Real Ron (10.0_dec);
    Real Cl = 0.0001_dec;
    Real Rl (1000.0_dec);

    /// Introduces the dynamics parameters
    Vector<Real> parameters(5);
    parameters[0] = amplitude; /// Amplitude of the input voltage, Vi
    parameters[1] = frequency; /// Sinusoid frequency, f
    parameters[2] = Ron; /// Diode resistance when on, Ron
    parameters[3] = Cl; /// Load capacitance, Cl
    parameters[4] = Rl; /// Load resistance, Rl

    RealConstant pi_c("pi",pi);

    /// Introduces the global parameters
    Real TIME_LIMIT = 1/frequency;
    //float TIME_LIMIT = 0.0042;
    Int TRAN_LIMIT = 1;
    double MAX_ENCL_RADIUS = 1.0;
    double MAX_STEP_SIZE = 1e-2/frequency.get_d();
    //float LOCK_TOGRID_TIME = 2.0/frequency;
    //double LOCK_TOGRID_TIME = 0.25/frequency;
    //Int MAX_GRID_DEPTH = 7;
    Bool ENABLE_SUBDIV=false;

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridAutomaton rectifier_automaton("rectifier_automaton");

    // Create the coordinates
    RealVariable t("t"); // Time
    RealVariable vi("Vin"); // Input voltage
    RealVariable vo("Vout"); // Output voltage

    Real one=1;

    /// Create the discrete states
    StringVariable rectifier("rectifier");
    DiscreteLocation offoff(rectifier|"offoff");
    DiscreteLocation onoff(rectifier|"onoff");
    DiscreteLocation offon(rectifier|"offon");
    DiscreteLocation onon(rectifier|"onon");

    /// Create the discrete events
    DiscreteEvent resettime(1);
    DiscreteEvent jump1(2), jump2(3), jump3(4);

    /// Create the resets

    /// Reset the time (t^=0,vi^=vi,vo^=vo)
    PrimedRealAssignments resettime_r( next({t,vi,vo}) = {Real(0),vi,vo} );
    /// Do nothing (t^=t,vi^=vi,vo^=vo)
    PrimedRealAssignments noop_r( next({t,vi,vo}) = {t,vi,vo} );

    /// Create the guards
    Real f=parameters[1];
    /// Guard for the reset of time (t>=1/f)
    ContinuousPredicate resettime_g( t>=1/f );
    /// Guard for the jump from onoff to offoff (vi-vo<=0)
    ContinuousPredicate onoff_offoff_g( vi<=vo );
    /// Guard for the jump from offon to offoff (-vi-vo<=0)
    ContinuousPredicate offon_offoff_g( vi+vo>=0 );
    /// Guard for the jump from offoff to onoff (vi-vo>=0)
    ContinuousPredicate offoff_onoff_g( vi>=vo );
    /// Guard for the jump from onon to onoff (-vi-vo<=0)
    ContinuousPredicate onon_onoff_g( vi+vo>=0 );
    /// Guard for the jump from offoff to offon (-vi-vo>=0)
    ContinuousPredicate offoff_offon_g( vi+vo<=0 );
    /// Guard for the jump from onon to offon (vi-vo<=0)
    ContinuousPredicate onon_offon_g( vi<=vo );
    /// Guard for the jump from offon to onon (vi-vo>=0)
    ContinuousPredicate offon_onon_g( vi>=vo );
    /// Guard for the jump from onoff to onon (-vi-vo>=0)
    ContinuousPredicate onoff_onon_g( vi+vo<=0 );

    /// Build the automaton

    /// Create the dynamics

    /// Dynamics for the case of both diodes being off
    /// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)
    DottedRealAssignments offoff_d( dot({t,vi,vo}) = {one,amplitude*2*pi*frequency*cos(2*pi*frequency*t),-vo/(Rl*Cl)} );
    /// Dynamics for the case of the first diode being on, the second being off
    /// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)+(vi-vo)/(Ron*Cl)
    RealExpressions onoff_d({one,amplitude*2*pi*frequency*cos(2*pi*frequency*t),-vo/(Rl*Cl)+(vi-vo)/(Ron*Cl)});
    /// Dynamics for the case of the first diode being off, the second being on
    /// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)-(vi+vo)/(Ron*Cl)
    RealExpressions offon_d({one,amplitude*2*pi*frequency*cos(2*pi*frequency*t),-vo/(Rl*Cl)+(vo-vi)/(Ron*Cl)});
    /// Dynamics for the case of both diodes being on
    /// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)-2*vo/(Ron*Cl)
    RealExpressions onon_d({one,amplitude*2*pi*frequency*cos(2*pi*frequency*t),-vo/(Rl*Cl)-2*vo/(Ron*Cl)});

    List<RealVariable> space( {t,vi,vo} );
    /// Locations
    rectifier_automaton.new_mode(offoff,offoff_d);
    rectifier_automaton.new_mode(onoff,dot(space)=onoff_d);
    rectifier_automaton.new_mode(offon,dot(space)=offon_d);
    rectifier_automaton.new_mode(onon,dot(space)=onon_d);
    /// OffOff events
    rectifier_automaton.new_transition(offoff,resettime,offoff,resettime_r,resettime_g,EventKind::URGENT);
    rectifier_automaton.new_transition(offoff,jump1,onoff,noop_r,offoff_onoff_g,EventKind::URGENT);
    rectifier_automaton.new_transition(offoff,jump2,offon,noop_r,offoff_offon_g,EventKind::URGENT);
    /// OnOff events
    rectifier_automaton.new_transition(onoff,resettime,onoff,resettime_r,resettime_g,EventKind::URGENT);
    rectifier_automaton.new_transition(onoff,jump1,offoff,noop_r,onoff_offoff_g,EventKind::URGENT);
    rectifier_automaton.new_transition(onoff,jump3,onon,noop_r,onoff_onon_g,EventKind::URGENT);
    /// OffOn events
    rectifier_automaton.new_transition(offon,resettime,offon,resettime_r,resettime_g,EventKind::URGENT);
    rectifier_automaton.new_transition(offon,jump1,offoff,noop_r,offon_offoff_g,EventKind::URGENT);
    rectifier_automaton.new_transition(offon,jump3,onon,noop_r,offon_onon_g,EventKind::URGENT);
    /// OnOn events
    rectifier_automaton.new_transition(onon,resettime,onon,resettime_r,resettime_g,EventKind::URGENT);
    rectifier_automaton.new_transition(onon,jump2,onoff,noop_r,onon_onoff_g,EventKind::URGENT);
    rectifier_automaton.new_transition(onon,jump3,offon,noop_r,onon_offon_g,EventKind::URGENT);


    /// Finished building the automaton

    cout << "Automaton = " << rectifier << endl << endl;

    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(rectifier_automaton);

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(MAX_ENCL_RADIUS);
    evolver.configuration().set_maximum_step_size(MAX_STEP_SIZE);
    evolver.configuration().set_enable_subdivisions(ENABLE_SUBDIV);
    std::cout <<  evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::OrbitType OrbitType;

    std::cout << "Computing evolution..." << std::endl;

    RealVariablesBox initial_box({t==0, vi==0, vo==Real(0.8_dec)*parameters[0]});
    HybridSet initial_set(offoff,initial_box);

//    ExactBoxType initial_box(3, 0.002836,0.002836, 3.110529,3.110529, 3.110529,3.110529);
//    HybridEnclosureType initial_enclosure(onoff,initial_box);

    std::cout << "Initial set=" << initial_set << std::endl;

    HybridTime evolution_time(TIME_LIMIT,TRAN_LIMIT);


    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final size="<<orbit.final().size()<<std::endl;

    Axes2d graphic_axes(0.0<=t<=1.0/parameters[1].get_d(),-parameters[0]<=vi<=parameters[0]);
    Axes2d graphic_axes2(-parameters[0]<=t<=parameters[0],2<=vi<=parameters[0]);

    std::cout << "Plotting results..." << std::flush;

    plot("rectifier_orbit_t_vin", graphic_axes, Colour(0.0,0.5,1.0), orbit);
    plot("rectifier_orbit_t_vout", graphic_axes, Colour(0.0,0.5,1.0), orbit);
    plot("rectifier_orbit_vin_vout", graphic_axes2, Colour(0.0,0.5,1.0), orbit);


/*
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = LOCK_TOGRID_TIME;
    analyser.parameters().maximum_grid_depth= MAX_GRID_DEPTH;
    rectifier_automaton.set_grid(Grid(Vector<FloatDP>({3, 0.25/dp[1], 1.0, 0.5})));
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[offoff]=initial_box;
    HybridTime reach_time(TIME_LIMIT,TRAN_LIMIT);

    std::cout << "Computing upper reach set... " << std::endl << std::flush;
    HybridGridTreeSet reach = analyser.upper_reach(rectifier,initial_set,reach_time);
    std::cout << "done." << std::endl;

    plot("rectifier_reach_t_vin", 0, 1, 3, graphic_box, Colour(0.0,0.5,1.0), reach, -1);
    plot("rectifier_reach_t_vout", 0, 2, 3, graphic_box, Colour(0.0,0.5,1.0), reach, -1);
    plot("rectifier_reach_vin_vout", 1, 2, 3, graphic_box2, Colour(0.0,0.5,1.0), reach, -1);
*/
}
