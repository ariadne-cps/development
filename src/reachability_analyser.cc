/*****************************************************************************************
 *            reachability_analyser.cc
 *
 *  Copyright  2006-10  Alberto Casagrande, Pieter Collins, Davide Bresolin, Luca Geretti
 *
 *****************************************************************************************/

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

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "exceptions.h"

#include "numeric.h"

#include "vector.h"
#include "matrix.h"

#include "box.h"
#include "list_set.h"
#include "grid_set.h"

#include "orbit.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"

#include "evolution_parameters.h"
#include "evolution_statistics.h"
#include "evolver_interface.h"

#include "discretiser.h"
#include "reachability_analyser.h"
#include "logging.h"

#include "graphics.h"


namespace Ariadne {

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser)
    : _parameters(new EvolutionParametersType())
    , _discretiser(discretiser.clone())
{
}






// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridReachabilityAnalyser::_upper_reach(const HybridAutomaton& sys,
                                         const HybridGridTreeSet& set,
                                         const HybridTime& time,
                                         const int accuracy) const
{
    Gr grid=sys.grid();
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        EnclosureType enclosure=this->_discretiser->enclosure(*iter);
        result.adjoin(this->_discretiser->reach(sys,enclosure,time,accuracy,UPPER_SEMANTICS));
    }
    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridAutomaton& sys,
                                          const HybridGridTreeSet& set,
                                          const HybridTime& time,
                                          const int accuracy) const
{
    Gr grid=sys.grid();
    GTS result(grid); GTS cells=set; cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        ARIADNE_LOG(6,"\t\t\t\t\tEvolving cell = "<<*iter<<"\n");
        EnclosureType enclosure=this->_discretiser->enclosure(*iter);
        result.adjoin(this->_discretiser->evolve(sys,enclosure,time,accuracy,UPPER_SEMANTICS));
    }
    ARIADNE_LOG(5,"\t\t\t\t_upper_evolve result size = "<<result.size()<<"\n");
    return result;
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomaton& sys,
                                                const HybridGridTreeSet& set,
                                                const HybridTime& time,
                                                const int accuracy) const
{
    ARIADNE_LOG(5,"\t\t\t\tHybridReachabilityAnalyser::_upper_reach_evolve(...)\n");
    Gr grid=sys.grid();
    std::pair<GTS,GTS> result=make_pair(GTS(grid),GTS(grid));
    GTS& reach=result.first; GTS& evolve=result.second;
    GTS cells=set; cells.mince(accuracy);

    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        EnclosureType enclosure=this->_discretiser->enclosure(*iter);
        GTS cell_reach, cell_final;
        make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(sys,enclosure,time,accuracy,UPPER_SEMANTICS);
        ARIADNE_LOG(7,"\t\t\t\t\t\tevolution reach size= "<<cell_reach.size()<<"\n");
        ARIADNE_LOG(7,"\t\t\t\t\t\tevolution final size= "<<cell_final.size()<<"\n");
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }
    ARIADNE_LOG(6,"\t\t\t\t\tfinal reach size = "<<reach.size()<<"\n");
    ARIADNE_LOG(6,"\t\t\t\t\tfinal evolve size = "<<evolve.size()<<"\n");
    ARIADNE_LOG(5,"Done.\n");
    return result;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::lower_evolve(...)\n");

	this->_discretiser->reset_lower_statistics(); // Resets the continuous statistics

    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    Gr grid(system.grid());
    GTS initial; GTS final;

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth+4));
        Vector<Float> origin = grid[loc_iter->first].origin();
        Box bbox = loc_iter->second.bounding_box();
        if (radius(bbox) > cell_radius) {
            // if bigger, map to the grid
            // First of all, test if the bounding box lies on cell boundaries or not
            for(uint i = 0 ; i < bbox.dimension() ; ++i) {
                // test if the i-th dimension is a singleton interval AND
                // if it lies on a cell boundary
                cell_radius = cell[i]/(1 << (grid_depth+4));
                Float intpart, fractpart;
                fractpart = modf((bbox[i].lower()-origin[i])/cell_radius,&intpart);
                if(bbox[i].singleton() && (fractpart == 0.0)) {
                    // the set lies on the boundary, shift the grid center by half cell size
                    origin[i]+=cell_radius/2.0;
                    grid[loc_iter->first].set_origin(origin);
                }
            }        
            // if bigger, map to the grid
            ARIADNE_LOG(6,"\t\t\t\t\tAdjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
            initial[loc_iter->first].adjoin_lower_approximation(loc_iter->second,grid_height,grid_depth+4);
        } else {
            ARIADNE_LOG(6,"\t\t\t\t\tComputing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_final=this->_discretiser->evolve(system,initial_enclosure,time,grid_depth,LOWER_SEMANTICS);
            final.adjoin(cell_final);            
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
        }
    }

    ARIADNE_LOG(5,"\t\t\t\tgrid="<<grid<<"\n");

    ARIADNE_LOG(5,"\t\t\t\tinitial.size()="<<initial.size()<<"\n");
    if(!initial.empty()) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing lower evolution from the grid.");
        for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
            ARIADNE_LOG(5,".");
            EnclosureType enclosure=this->_discretiser->enclosure(*bs_iter);
            GTS cell_final=this->_discretiser->evolve(system,enclosure,time,grid_depth,LOWER_SEMANTICS);
            final.adjoin(cell_final);
        }
    }
    ARIADNE_LOG(4,"\n");

	// Copies the largest evolution time and steps to the statistics
	this->_statistics->lower().largest_evol_time = this->_discretiser->statistics().lower().largest_evol_time;
	this->_statistics->lower().largest_evol_steps = this->_discretiser->statistics().lower().largest_evol_steps;

    return final;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::lower_reach(...)\n");

	this->_statistics->lower().reset(); // Resets the discrete statistics
	this->_discretiser->reset_lower_statistics(); // Resets the continuous statistics

    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    Gr grid(system.grid());
    GTS initial; GTS reach;

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth+4));
        Vector<Float> origin = grid[loc_iter->first].origin();
        Box bbox = loc_iter->second.bounding_box();
        if (radius(bbox) > cell_radius) {
            // if bigger, map to the grid
            // First of all, test if the bounding box lies on cell boundaries or not
            for(uint i = 0 ; i < bbox.dimension() ; ++i) {
                // test if the i-th dimension is a singleton interval AND
                // if it lies on a cell boundary
                cell_radius = cell[i]/(1 << (grid_depth+4));
                Float intpart, fractpart;
                fractpart = modf((bbox[i].lower()-origin[i])/cell_radius,&intpart);
                if(bbox[i].singleton() && (fractpart == 0.0)) {
                    // the set lies on the boundary, shift the grid center by half cell size
                    origin[i]+=cell_radius/2.0;
                    grid[loc_iter->first].set_origin(origin);
                }
            }        
            // if bigger, map to the grid
            ARIADNE_LOG(6,"\t\t\t\t\tAdjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
            initial[loc_iter->first].adjoin_lower_approximation(loc_iter->second,grid_height,grid_depth+4);
        } else {
            ARIADNE_LOG(6,"\t\t\t\t\tComputing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach=this->_discretiser->reach(system,initial_enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);            
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
        }
    }

    ARIADNE_LOG(5,"\t\t\t\tgrid="<<grid<<"\n");

    ARIADNE_LOG(5,"\t\t\t\tinitial.size()="<<initial.size()<<"\n");
    if(!initial.empty()) {
        ARIADNE_LOG(5,"\t\t\t\tComputing lower reach set from the grid...");
        for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
            ARIADNE_LOG(5,".");
            EnclosureType enclosure=this->_discretiser->enclosure(*bs_iter);
            GTS cell_reach = this->_discretiser->reach(system,enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);
        }
    }

	// Copies the reached region
	this->_statistics->lower().reach = reach;
	// Copies the largest evolution time and steps to the statistics
	this->_statistics->lower().largest_evol_time = this->_discretiser->statistics().lower().largest_evol_time;
	this->_statistics->lower().largest_evol_steps = this->_discretiser->statistics().lower().largest_evol_steps;

	return reach;
}



std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::lower_reach_evolve(...)\n");

	this->_statistics->lower().reset(); // Resets the discrete statistics
	this->_discretiser->reset_lower_statistics(); // Resets the continuous statistics

    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;

    Gr grid=system.grid();
    GTS initial;
    GTS reach; GTS evolve;

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth+4));
        Vector<Float> origin = grid[loc_iter->first].origin();
        Box bbox = loc_iter->second.bounding_box();
        if (radius(bbox) > cell_radius) {
            // if bigger, map to the grid
            // First of all, test if the bounding box lies on cell boundaries or not
            for(uint i = 0 ; i < bbox.dimension() ; ++i) {
                // test if the i-th dimension is a singleton interval AND
                // if it lies on a cell boundary
                cell_radius = cell[i]/(1 << (grid_depth+4));
                Float intpart, fractpart;
                fractpart = modf((bbox[i].lower()-origin[i])/cell_radius,&intpart);
                if(bbox[i].singleton() && (fractpart == 0.0)) {
                    // the set lies on the boundary, shift the grid center by half cell size
                    origin[i]+=cell_radius/2.0;
                    grid[loc_iter->first].set_origin(origin);
                }
            }        
            // if bigger, map to the grid
            ARIADNE_LOG(6,"\t\t\t\t\tAdjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));
            initial[loc_iter->first].adjoin_lower_approximation(loc_iter->second,grid_height,grid_depth+4);
        } else {
            ARIADNE_LOG(6,"\t\t\t\t\tComputing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach, cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
            initial.insert(make_pair(loc_iter->first,grid[loc_iter->first]));           
        }
    }

    ARIADNE_LOG(5,"\t\t\t\tgrid="<<grid<<"\n");

    ARIADNE_LOG(5,"\t\t\t\tinitial.size()="<<initial.size()<<"\n");
    if(!initial.empty()) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing lower evolution from the grid.");
        for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
            ARIADNE_LOG(5,".");
            EnclosureType enclosure=this->_discretiser->enclosure(*bs_iter);
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final) = this->_discretiser->evolution(system,enclosure,time,grid_depth,LOWER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }

	// Copies the reached region
	this->_statistics->lower().reach = reach;
	// Copies the largest evolution time and steps to the statistics
	this->_statistics->lower().largest_evol_time = this->_discretiser->statistics().lower().largest_evol_time;
	this->_statistics->lower().largest_evol_steps = this->_discretiser->statistics().lower().largest_evol_steps;

    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::upper_evolve(...)\n");
 
	this->_statistics->upper().reset(); // Resets the discrete statistics
	this->_discretiser->reset_upper_statistics(); // Resets the continuous statistics

    Gr grid=system.grid();
    GTS evolve(grid), initial(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remaining_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remaining_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remaining_time(remaining_time,discrete_steps);
    ARIADNE_LOG(5,"\t\t\t\treal_time="<<real_time<<"\n");
    ARIADNE_LOG(5,"\t\t\t\ttime_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(5,"\t\t\t\tcomputing first reachability step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(5,"\t\t\t\tAdjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,grid_depth);
        } else {
            ARIADNE_LOG(5,"\t\t\t\tComputing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_final=this->_discretiser->evolve(system,initial_enclosure,time,grid_depth,UPPER_SEMANTICS);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing evolution from the grid...\n");
        ARIADNE_LOG(6,"\t\t\t\t\tinitial_evolve.size()="<<initial.size()<<"\n");
        initial=this->_upper_evolve(system,initial,hybrid_lock_to_grid_time,grid_depth);
        evolve.adjoin(initial);
    }

	// Adds the largest evolution time and steps to the statistics, then resets such statistics
	this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time;
	this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps;
	this->_discretiser->reset_upper_largest_evol_statistics();

    // time steps evolution loop
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing "<<i+1<<"-th reachability step...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
		this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,evolve.size()); // Updates the largest intermediate size
		this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time; // Adds the largest evolution time to the statistics
		this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps; // Adds the largest evolution steps to the statistics
		this->_discretiser->reset_upper_largest_evol_statistics(); // Resets the continuous statistics related to the largest evolution time/steps
    }

	this->_statistics->upper().total_locks = time_steps+1; // Sets the number of locks (time_steps + the "initial" lock)

    ARIADNE_LOG(5,"\t\t\t\tremaining_time="<<remaining_time<<"\n");
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing evolution for remaining time...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_remaining_time,grid_depth);
		this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,evolve.size()); // Updates the largest intermediate size
		this->_statistics->upper().total_locks++; // Increases the total locks counter
		this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time; // Adds the largest evolution time to the statistics
		this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps; // Adds the largest evolution steps to the statistics
		this->_discretiser->reset_upper_largest_evol_statistics(); // Resets the continuous statistics related to the largest evolution time/steps
    }
    evolve.recombine();
    ARIADNE_LOG(6,"\t\t\t\t\tfinal_evolve.size()="<<evolve.size()<<"\n");
    return evolve;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::upper_reach(system,set,time)\n");

	GTS reach, evolve;
	make_lpair(reach,evolve) = this->upper_reach_evolve(system, initial_set, time); // Runs the upper_reach_evolve routine on its behalf

    return reach; // Returns the reached region only
}



std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::upper_reach_evolve(system,set,time)\n");
    ARIADNE_LOG(5,"\t\t\t\tinitial_set="<<initial_set<<"\n");

	this->_statistics->upper().reset(); // Reset the discrete statistics
	this->_discretiser->reset_upper_statistics(); // Reset the continuous statistics

    Gr grid=system.grid();
    GTS found(grid),evolve(grid),reach(grid),initial(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remaining_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remaining_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remaining_time(remaining_time,discrete_steps);
    ARIADNE_LOG(5,"\t\t\t\treal_time="<<real_time<<"\n");
    ARIADNE_LOG(5,"\t\t\t\ttime_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(5,"\t\t\t\tcomputing first reachability step...\n");
    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (grid_depth));
        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(6,"\t\t\t\t\tAdjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,grid_depth);
        } else {
            ARIADNE_LOG(6,"\t\t\t\t\tComputing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,time,grid_depth,UPPER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }
    if(!initial.empty()) {
        ARIADNE_LOG(5,"\t\t\t\t\tcomputing evolution from the grid...\n");
        ARIADNE_LOG(6,"\t\t\t\t\tinitial_evolve.size()="<<initial.size()<<"\n");
        make_lpair(found,initial)=this->_upper_reach_evolve(system,initial,hybrid_lock_to_grid_time,grid_depth);
        evolve.adjoin(initial);
		this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,initial.size()); // Updates the largest intermediate size
        reach.adjoin(found);
    }

	// Adds the largest evolution time and steps to the statistics, then resets such statistics
	this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time;
	this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps;
	this->_discretiser->reset_upper_largest_evol_statistics();

    // time steps evolution loop        
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(6,"\t\t\t\t\tfound.size()="<<found.size()<<"\n");
        ARIADNE_LOG(6,"\t\t\t\t\tevolve.size()="<<evolve.size()<<"\n");
        reach.adjoin(found);
		this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,evolve.size()); // Updates the largest intermediate size
		this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time; // Adds the largest evolution time to the statistics
		this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps; // Adds the largest evolution steps to the statistics
		this->_discretiser->reset_upper_largest_evol_statistics(); // Resets the continuous statistics related to the largest evolution time/steps
        ARIADNE_LOG(5,"\t\t\t\tfound "<<found.size()<<" cells.\n");
    }

	this->_statistics->upper().total_locks = time_steps+1; // Sets the number of locks (time_steps + the "initial" lock)

    ARIADNE_LOG(5,"\t\t\t\tremaining_time="<<remaining_time<<"\n");
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(5,"\t\t\t\tcomputing evolution for the remaining time...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_remaining_time,grid_depth);
        reach.adjoin(found);
		this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,evolve.size()); // Updates the largest intermediate size
		this->_statistics->upper().total_locks++; // Increases the total locks counter
		this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time; // Adds the largest evolution time to the statistics
		this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps; // Adds the largest evolution steps to the statistics
		this->_discretiser->reset_upper_largest_evol_statistics(); // Resets the continuous statistics related to the largest evolution time/steps
    }

    reach.recombine();
	this->_statistics->upper().reach = reach;
    ARIADNE_LOG(5,"\t\t\t\treach="<<reach<<"\n");
    evolve.recombine();
    ARIADNE_LOG(5,"\t\t\t\tevolve="<<evolve<<"\n");
    return std::make_pair(reach,evolve);
}




HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set) const
{
    ARIADNE_LOG(4,"\t\t\tHybridReachabilityAnalyser::chain_reach(system,initial_set)\n");

	this->_statistics->upper().reset(); // Resets the discrete statistics
	this->_discretiser->reset_upper_statistics(); // Resets the continuous statistics

    HybridBoxes bounding_domain = this->_parameters->bounding_domain;
    Float transient_time = this->_parameters->transient_time;
    int transient_steps = this->_parameters->transient_steps;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    int lock_to_grid_steps = this->_parameters->lock_to_grid_steps;
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;

    ARIADNE_LOG(6,"\t\t\t\t\ttransient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(6,"\t\t\t\t\tlock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");
    ARIADNE_LOG(6,"\t\t\t\t\tbounding_domain="<<bounding_domain<<"\n");
    ARIADNE_LOG(6,"\t\t\t\t\tinitial_set="<<initial_set<<"\n");

	// Checks consistency of the bounding domain in respect to the state space
	HybridSpace hspace = system.state_space();
	// If the DiscreteState was not found or otherwise if the continuous space sizes mismatch, throws an error
	for (HybridSpace::locations_const_iterator hs_it = hspace.locations_begin(); hs_it != hspace.locations_end(); ++hs_it) 
	{
		if (bounding_domain.find(hs_it->first) == bounding_domain.end())
		{
			ARIADNE_FAIL_MSG("Error: the system state space and the bounding domain space do not match on the discrete space.");		
		}		
		else if (hs_it->second != bounding_domain[hs_it->first].size())
			ARIADNE_FAIL_MSG("Error: the system state space and the bounding domain space do not match on the continuous space.");		
	}

    Gr grid=system.grid();
    GTS bounding(grid), evolve(grid), reach(grid), initial(grid), found(grid);

    bounding.adjoin_outer_approximation(bounding_domain,maximum_grid_depth); 
	bounding.recombine();
	HybridBoxes bounding_box = bounding.bounding_box(); // Used for the restriction check
    ARIADNE_LOG(6,"\t\t\t\t\tbounding_size="<<bounding.size()<<"\n");

    if(transient_time <= 0.0 || transient_steps <= 0) {
        transient_time = lock_to_grid_time;
        transient_steps = lock_to_grid_steps;
    }
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    ARIADNE_LOG(5,"\t\t\t\tComputing first evolution step...\n");

    // For each location, test if the radius of the set is smaller than the grid cell
    for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
        loc_iter!=initial_set.locations_end(); ++loc_iter) 
    {
        Vector<Float> cell = grid[loc_iter->first].lengths();
        Float cell_radius = (min(cell))/(1 << (maximum_grid_depth));

        if (radius(loc_iter->second.bounding_box()) > cell_radius) {
            // if bigger, map to the grid
            ARIADNE_LOG(6,"\t\t\t\t\tAdjoining initial set for location "<<loc_iter->first<<" to the grid...\n");
            initial[loc_iter->first].adjoin_outer_approximation(loc_iter->second,maximum_grid_depth);
        } else {
            ARIADNE_LOG(6,"\t\t\t\t\tComputing evolution for initial set in location "<<loc_iter->first<<" directly...\n");
            // if smaller, compute the evolution directly
            EnclosureType initial_enclosure(loc_iter->first,ContinuousEnclosureType(loc_iter->second));
            GTS cell_reach,cell_final;
            make_lpair(cell_reach,cell_final)=this->_discretiser->evolution(system,initial_enclosure,hybrid_transient_time,maximum_grid_depth,UPPER_SEMANTICS);
            reach.adjoin(cell_reach);
            evolve.adjoin(cell_final);
        }
    }

    if(!initial.empty()) {
        ARIADNE_LOG(5,"\t\t\t\tComputing evolution on the grid...\n")
        make_lpair(found,initial)=this->_upper_reach_evolve(system,initial,hybrid_transient_time,maximum_grid_depth);
        ARIADNE_LOG(6,"\t\t\t\t\tfound "<<found.size()<<" cells.\n");
        reach.adjoin(found);
        evolve.adjoin(initial);
    }

	// Adds the largest evolution time and steps to the statistics, then resets such statistics
	this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time;
	this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps;
	this->_discretiser->reset_upper_largest_evol_statistics();
	 
    // If the evolve region is not a subset of the bounding region, the region will be restricted (NOTE: for efficiency, only performed if the region is currently considered unrestricted)
	if (!this->_statistics->upper().has_restriction_occurred)
		if (!evolve.subset(bounding_box)) this->_statistics->upper().has_restriction_occurred = true;

    evolve.restrict(bounding);
	this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,evolve.size()); // Updates the largest intermediate size

    ARIADNE_LOG(5,"\t\t\t\tfound "<<reach.size()<<" cells, of which "<<evolve.size()<<" are new.\n");   
    
    ARIADNE_LOG(5,"\t\t\t\tComputing recurrent evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

	// While the final set has new cells in respect to the previous reach set, process them and increase the number of locks (starting from 1 due to the initial transient phase)
    for (this->_statistics->upper().total_locks = 1; !evolve.empty(); this->_statistics->upper().total_locks++) {

        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(6,"\t\t\t\t\tfound.size()="<<found.size()<<"\n");
        ARIADNE_LOG(6,"\t\t\t\t\tevolve.size()="<<evolve.size()<<"\n");

		// If the evolve region is not a subset of the bounding region, the region will be restricted (NOTE: for efficiency, only performed if the region is currently considered unrestricted)
		if (!this->_statistics->upper().has_restriction_occurred)
			if (!evolve.subset(bounding_box)) this->_statistics->upper().has_restriction_occurred = true;

        evolve.remove(reach);
        evolve.restrict(bounding);
        reach.adjoin(found);

		this->_statistics->upper().largest_intermediate_size = max(this->_statistics->upper().largest_intermediate_size,evolve.size()); // Update the largest intermediate size
		this->_statistics->upper().largest_evol_time += this->_discretiser->statistics().upper().largest_evol_time; // Adds the largest evolution time to the statistics
		this->_statistics->upper().largest_evol_steps += this->_discretiser->statistics().upper().largest_evol_steps; // Adds the largest evolution steps to the statistics
		this->_discretiser->reset_upper_largest_evol_statistics(); // Resets the continuous statistics related to the largest evolution time/steps

        ARIADNE_LOG(6,"\t\t\t\t\tfound "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
    }
    reach.recombine();

    // If the evolve region is not a subset of the bounding region, the region will be restricted (NOTE: for efficiency, only performed if the region is currently considered unrestricted)
	if (!this->_statistics->upper().has_restriction_occurred)
		if (!evolve.subset(bounding_box)) this->_statistics->upper().has_restriction_occurred = true;

    reach.restrict(bounding);
	this->_statistics->upper().reach = reach;

    return reach;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const HybridBoxes& bounding_set) const
{
	// Assigns the input bounding_set to the bounding domain
	this->_parameters->bounding_domain = bounding_set;

	// Returns the result of chain_reach with implicit bounding set
	return this->chain_reach(system,initial_set);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
viable(const SystemType& system,
       const HybridImageSet& bounding_set) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


bool 
HybridReachabilityAnalyser::
_safe(const SystemType& system, 
	  const HybridImageSet& initial_set, 
	  const HybridBoxes& safe_box)
{
	ARIADNE_LOG(4,"\t\t\tSafety analysis...\n");

	HybridGridTreeSet reach = chain_reach(system,initial_set); // Perform the chain reachability analysis	

	// If the reached region was not restricted and the maximum enclosure bounds have not been reached while not subdividing, a verification is feasible
	if (!this->_statistics->upper().has_restriction_occurred && !(this->_discretiser->statistics().upper().has_max_enclosure_been_reached && !this->_discretiser->parameters().enable_subdivisions))
	{
		// If the reached region is definitely inside the hybrid safe box, the result is safe 
		bool result = definitely(reach.subset(safe_box));

		ARIADNE_LOG(4, (result ? "\t\t\tSafe.\n" : "\t\t\tNot safe.\n") );

		return result;
	}
	// Otherwise notify and return false
	else
	{
		ARIADNE_LOG(4,"\t\t\tNot checked due to: ");
		if (this->_statistics->upper().has_restriction_occurred)
			ARIADNE_LOG(4,"<domain bounds> ");
		if (this->_discretiser->statistics().upper().has_max_enclosure_been_reached && !this->_discretiser->parameters().enable_subdivisions)
			ARIADNE_LOG(4,"<enclosure bounds>");
		ARIADNE_LOG(4,"\n");
		return false;
	}
}


bool 
HybridReachabilityAnalyser::
_unsafe(const SystemType& system, 
		const HybridImageSet& initial_set, 
		const HybridBoxes& safe_box)
{
	ARIADNE_LOG(4,"\t\t\tUnsafety analysis...\n");

	// Get the size of the continuous space (NOTE: assumed equal for all locations)
	const uint css = system.state_space().locations_begin()->second; 
	// Create the evolution time from the upper approximation obtained in the upper case
	HybridTime lrt = HybridTime(this->_statistics->upper().largest_evol_time,this->_statistics->upper().largest_evol_steps); 

	// Get the sizes of the upper reached set for each location

		// Get a copy of the upper reached region
		HybridGridTreeSet upperreach = this->_statistics->upper().reach;
		// The sizes of the reached region
		std::map<DiscreteState,uint> sizes;
		// For each location, obtain the size of the reach set
		for(HybridGridTreeSet::locations_iterator loc_iter=upperreach.locations_begin();
			loc_iter!=upperreach.locations_end(); ++loc_iter) 
		{
			// Mince the region to the maximum depth
			loc_iter->second.mince(this->_parameters->maximum_grid_depth);
			// Store the (location,size) pair
			sizes.insert(make_pair<DiscreteState,uint>(loc_iter->first,loc_iter->second.size()));
		}
		// Assign the sizes
		this->_discretiser->parameters().hybrid_working_sets_size_limit = sizes;

	HybridGridTreeSet lowerreach = lower_reach(system,initial_set,lrt);

	// Create a safe set enlarged for the falsification		
	HybridBoxes safe_box_enl = safe_box;
	for (HybridBoxes::iterator hbx_it = safe_box_enl.begin(); hbx_it != safe_box_enl.end(); hbx_it++)
		hbx_it->second.widen();

	// Get a copy of the grid
	HybridGrid hg = system.grid();
	// Get the safe box enlarged by a grid cell and half a maximum enclosure cell
	for (HybridGrid::locations_const_iterator hg_it = hg.locations_begin(); hg_it != hg.locations_end(); hg_it++)
		for (uint i=0;i<css;i++)
			safe_box_enl[hg_it->first][i] += Interval(-this->_discretiser->statistics().lower().largest_enclosure_cell[i]/2,this->_discretiser->statistics().lower().largest_enclosure_cell[i]/2)+Interval(-hg_it->second.lengths()[i],hg_it->second.lengths()[i])/(1<<this->_parameters->maximum_grid_depth);

	// If the reach is definitely not included in the enlarged safe box, then it is unsafe
	bool result = definitely(!lowerreach.subset(safe_box_enl));

	ARIADNE_LOG(4, (result ? "\t\t\tUnsafe.\n" : "\t\t\tNot unsafe.\n") );

	return result;
}


tribool
HybridReachabilityAnalyser::
verify(const SystemType& system,
       const HybridImageSet& initial_set,
       const HybridBoxes& safe_box)
{
		ARIADNE_LOG(3, "\t\tVerification...\n");

		// Perform the safety analysis
		if (this->_safe(system,initial_set,safe_box)) 
		{
			ARIADNE_LOG(3, "\t\tSafe.\n");
			return true;
		}
	
		ARIADNE_LOG(4, "\t\t\tChain reach largest evolution time: " << this->_statistics->upper().largest_evol_time << "\n");
		ARIADNE_LOG(4, "\t\t\tChain reach largest evolution steps: " << this->_statistics->upper().largest_evol_steps << "\n");

		// Perform the unsafety analysis
		if (this->_unsafe(system,initial_set,safe_box)) 
		{
			ARIADNE_LOG(3, "\t\tUnsafe.\n");
			return false;
		}

		ARIADNE_LOG(4, "\t\t\tLower reach largest evolution time: " << this->_statistics->lower().largest_evol_time << "\n");
		ARIADNE_LOG(4, "\t\t\tLower reach largest evolution steps: " << this->_statistics->lower().largest_evol_steps << "\n");
		ARIADNE_LOG(4, "\t\t\tLower reach largest enclosure cell: " << this->_discretiser->statistics().lower().largest_enclosure_cell << "\n");

		ARIADNE_LOG(3, "\t\tIndeterminate.\n");

		// Return indeterminate if both failed
		return indeterminate;
}


HybridFloatVector 
HybridReachabilityAnalyser::
_getDomainHMAD(const HybridAutomaton& system) const
{ 
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;
	// The variable for the result, correctly initialized
	Vector<Float> mad = Vector<Float>(css);
    // The variable for the hybrid maximum absolute derivatives
	HybridFloatVector hmad;

    // For each mode
	for (set<DiscreteMode>::const_iterator it = system.modes().begin(); it != system.modes().end(); it++)
	{
		// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
		der = it->dynamic()(this->_parameters->bounding_domain.find(it->location())->second); 

		// Gets the maximum absolute derivatives
		for (uint i=0;i<css;i++)
			mad[i] = abs(der[i]).upper();

		// Inserts the values for the location
		hmad.insert(pair<DiscreteState,Vector<Float> >(it->location(),mad)); 
	}

	// Returns
	return hmad;
}


HybridFloatVector 
HybridReachabilityAnalyser::
_getReachHMAD(const HybridAutomaton& system) const
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;
    // The variable for the result
	HybridFloatVector hmad;

	// Gets the dynamic for the given DiscreteState
	for (set<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++)
	{
		// Inserts the corresponding pair, initialized with zero maximum absolute derivatives
		hmad.insert(pair<DiscreteState,Vector<Float> >(modes_it->location(),Vector<Float>(css)));

		// For each of its hybrid cells			
		for (GridTreeSet::const_iterator cells_it = this->_statistics->upper().reach[modes_it->location()].begin(); 
										 cells_it != this->_statistics->upper().reach[modes_it->location()].end(); 
										 cells_it++)
		{
			// Gets the derivative bounds
			der = modes_it->dynamic()(cells_it->box());

			// For each variable, sets the maximum value
			for (uint i=0;i<css;i++)
				hmad[modes_it->location()][i] = max(hmad[modes_it->location()][i], abs(der[i]).upper()); 
		}			
	}

	// Returns
	return hmad;
}


void 
HybridReachabilityAnalyser::
_setMaximumStepSize(const HybridFloatVector& hmad, 
					const HybridGrid& hgrid)
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initializes the maximum step size
	Float mss = std::numeric_limits<double>::infinity();
	// For each couple DiscreteState,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) 
	{
		// For each dimension of the space
		for (uint i=0;i<css;i++)
			mss = min(mss,hgrid[hfv_it->first].lengths()[i]/(1<<this->_parameters->maximum_grid_depth)/hfv_it->second[i]);
	}

	// Assigns
	this->_discretiser->parameters().maximum_step_size = mss;
}

void
HybridReachabilityAnalyser::
_setMaximumEnclosureCell(const HybridGrid& hgrid)
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hgrid.locations_begin()->second.lengths().size();

	// Initializes the result
	Vector<Float> mec(css);

	// For each location and dimension of the space
	for (HybridGrid::locations_const_iterator hg_it = hgrid.locations_begin(); hg_it != hgrid.locations_end(); hg_it++)
		for (uint i=0;i<css;i++)
			if (hg_it->second.lengths()[i] > mec[i])
				mec[i] = hg_it->second.lengths()[i];

	// Scales the cell in respect to the maximum grid depth
	for (uint i=0;i<css;i++)
		mec[i] /= (1<<this->_parameters->maximum_grid_depth);

	// Assigns
	this->_discretiser->parameters().maximum_enclosure_cell = 2*mec;
}


HybridGrid 
HybridReachabilityAnalyser::
_getHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid grid
	HybridGrid hg;	

	// For each couple DiscreteState,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Initialize the maximum derivative/domain_width ratio
		Float maxratio = 0.0;
		// Initialize the gridlengths
		Vector<Float> gridlengths(css);
		// Initialize the minimum length to use for zero-derivative variables (will be the minimum among both grid lengths and domain widths)
		Float minlength = std::numeric_limits<double>::infinity();

		// For each dimension of the continuous space
		for (uint i=0;i<css;i++)
		{	
			maxratio = max(maxratio,hfv_it->second[i]/this->_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Get the largest ratio
			minlength = min(minlength,this->_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Get the smallest domain width
		}

		// Assign the lengths and check the minimum length to be assigned
		for (uint i=0;i<css;i++)
		{
			if (hfv_it->second[i] > 0) // If the derivative is greater than zero
			{
				gridlengths[i] = hfv_it->second[i]/maxratio; // Assign the length
				minlength = min(minlength,gridlengths[i]); // Reduce the minimum length by comparing it to the current grid length
			}
		}

		// Assign the minimum length for zero-derivative variables
		for (uint i=0;i<css;i++)
			if (gridlengths[i] == 0) // If it has zero grid length, assign minlength
				gridlengths[i] = minlength;

		// Assign the grid, centered on the origin
		hg[hfv_it->first] = Grid(Vector<Float>(css),gridlengths);

	}

	// Return
	return hg;
}


HybridGrid 
HybridReachabilityAnalyser::
_getEqualizedHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the maximum derivative/domain_width ratio
	Float maxratio = 0.0;
	// Initialize the minimum length to use for zero-derivative variables
	Float minlength = std::numeric_limits<double>::infinity();
	// Initialize the gridlengths
	Vector<Float> gridlengths(css);		
	// Initialize the maximum absolute derivatives
	Vector<Float> mad(css);

	// Get the maximum absolute derivatives
	for (uint i=0;i<css;i++)
		// For each couple DiscreteState,Vector<Float>
		for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
			mad[i] = max(mad[i],hfv_it->second[i]);

	// For each dimension
	for (uint i=0;i<css;i++)
		// For each couple DiscreteState,Vector<Float>
		for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
		{
			maxratio = max(maxratio,mad[i]/this->_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Check for maximum ratio
			minlength = min(minlength,this->_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Get the minimum domain length
		}

	// Assign the lengths and checks the minimum length to be assigned
	for (uint i=0;i<css;i++)
	{
		// If the derivative is greater than zero
		if (mad[i] > 0)
		{	
			gridlengths[i] = mad[i]/maxratio; // Assign the length
			minlength = min(minlength,gridlengths[i]); // Get the minimum between the current minimum and the current length
		}		
	}

	// Correct the lengths that are zero
	for (uint i=0;i<css;i++)
		// If it has zero grid length, assign minlength
		if (gridlengths[i] == 0)
			gridlengths[i] = minlength;
	
	// Save the grid, centered on the origin
	Grid gr(Vector<Float>(css),gridlengths);

	// Initialize the hybrid grid
	HybridGrid hg;
	// Populate it
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
		hg[hfv_it->first] = gr;		
	
	// Return
	return hg;
} 


void 
HybridReachabilityAnalyser::
_setInitialParameters(SystemType& system, const HybridBoxes& domain)
{
	// Fixed parameters
	this->_discretiser->parameters().enable_subdivisions = true;
	this->_discretiser->parameters().enable_set_model_reduction = true;
	this->_discretiser->parameters().enable_working_sets_pruning = true;
	this->_parameters->transient_time = 1e10;
	this->_parameters->transient_steps = 1;
	this->_parameters->lock_to_grid_time = 1e10;		
	this->_parameters->lock_to_grid_steps = 1;
	this->_parameters->bounding_domain = domain; // Set the domain (IMPORTANT: must be done before using _getDomainHMAD and _getEqualizedHybridGrid)

	this->_parameters->maximum_grid_depth = 0; // Initial value, incremented at each iteration

	HybridFloatVector hmad = this->_getDomainHMAD(system); // Evaluate the maximum absolute derivatives from the domain

	system.set_grid(this->_getEqualizedHybridGrid(hmad)); // Initial grid
	this->_setMaximumStepSize(hmad,system.grid()); // Initial maximum step size
	this->_setMaximumEnclosureCell(system.grid()); // Initial maximum enclosure cell
}


void 
HybridReachabilityAnalyser::
_adaptParameters(SystemType& system)
{
		// Evaluate the maximum absolute derivatives
		HybridFloatVector hmad = (this->_statistics->upper().reach.size() > 0) ? this->_getReachHMAD(system) : this->_getDomainHMAD(system); 

		this->_parameters->maximum_grid_depth++; // Increase the maximum grid depth

		system.set_grid(this->_getEqualizedHybridGrid(hmad)); // New grid
		this->_setMaximumStepSize(hmad,system.grid()); // New maximum step size
		this->_setMaximumEnclosureCell(system.grid()); // New maximum enclosure cell
}


tribool 
HybridReachabilityAnalyser::
verify_iterative(SystemType& system, 
				 const HybridImageSet& initial_set, 
				 const HybridBoxes& safe_box, 
				 const HybridBoxes& domain)
{
	ARIADNE_LOG(2,"\n\tIterative verification...\n");

	// Limit size for the chain reach
	const uint MAX_REACH_SIZE = 150000;
	// Limit maximum number of working sets for the lower reach
	const uint MAX_WORKING_SETS_TOTAL = 50000;

	// Initialize the seed for internal random number generation
	srand(time(NULL));

	// Save the folder name as a function of the automaton name and of the current timestamp, then create the main folder and the verification run folder
	time_t mytime;
	time(&mytime);
	string foldername = system.name()+".png";
	mkdir(foldername.c_str(),0777);
	foldername = foldername+"/"+asctime(localtime(&mytime));
	mkdir(foldername.c_str(),0777);
	char mgd_char[10];
	string filename;

	// Set the initial parameters
	this->_setInitialParameters(system, domain);

    while(1)
	{ 
		/// Print some information on the current iteration
		sprintf(mgd_char,"%i",this->_parameters->maximum_grid_depth);
		ARIADNE_LOG(2, "\tDEPTH " << this->_parameters->maximum_grid_depth << "\n"); 
		ARIADNE_LOG(3, "\n\t\tMaximum step size: " << this->_discretiser->parameters().maximum_step_size << "\n");
		ARIADNE_LOG(3, "\t\tMaximum enclosure cell: " << this->_discretiser->parameters().maximum_enclosure_cell << "\n");

		// Perform the verification
		tribool result = this->verify(system,initial_set,safe_box);
		// Return the result, if it is not indeterminate
		if (!indeterminate(result)) return result;

		// Plot the reached region (if no definite result has been obtained)
		filename = "upper-";
		plot(foldername,filename + mgd_char, this->_statistics->upper().reach);
		// Plot the reached region (if no definite result has been obtained)
		filename = "lower-";
		plot(foldername,filename + mgd_char, this->_statistics->lower().reach); 

		// Check if the chain reach size is sufficiently low, and that the lower reach working sets total is sufficiently low, in order to stop the verification
		if (this->_statistics->upper().reach.size() > MAX_REACH_SIZE || this->_discretiser->statistics().lower().largest_working_sets_total > MAX_WORKING_SETS_TOTAL) break;
		
		// Adapt the parameters for the next iteration
		this->_adaptParameters(system);
    }

	// Return indeterminate
	return indeterminate;
}

} // namespace Ariadne
