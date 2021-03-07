/***************************************************************************
 *            reachability_analyser.tpl.hpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins, Davide Bresolin
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

#include "function/functional.hpp"

#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>


#include "utility/exceptions.hpp"

#include "numeric/numeric.hpp"

#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"

#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "geometry/grid_paving.hpp"

#include "solvers/integrator.hpp"
#include "solvers/solver.hpp"
#include "geometry/function_set.hpp"

#include "dynamics/reachability_analyser.hpp"

#include "output/logging.hpp"
#include "output/graphics.hpp"
#include "output/progress_indicator.hpp"
#include "solvers/linear_programming.hpp"

namespace Ariadne {

class IteratedMap;

OutputStream& operator<<(OutputStream& os, const ChainOverspillPolicy& policy);



template<class SYS>
ReachabilityAnalyser<SYS>::
ReachabilityAnalyser(const EvolverType& evolver)
    : _system(evolver.system().clone())
    , _evolver(evolver.clone())
    , _configuration(new ConfigurationType(*this))
{
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
clone() const
    -> ReachabilityAnalyser<SYS>*
{
    return new ReachabilityAnalyser(*this);
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
_reach_evolve_resume(const ListSet<EnclosureType>& initial_enclosures,
                     const TimeType& time,
                     const Nat accuracy,
                     ListSet<EnclosureType>& evolve_enclosures,
                     Semantics semantics,
                     const EvolverType& evolver) const
    -> Pair<StorageType,StorageType>
{
    ARIADNE_LOG_SCOPE_CREATE;
    const GridType& grid=this->_configuration->grid();
    Pair<StorageType,StorageType> result=make_pair(StorageType(grid,this->system()),StorageType(grid,this->system()));
    StorageType& reach_cells=result.first; StorageType& evolve_cells=result.second;

    for(auto initial_enclosure : initial_enclosures) {
        ListSet<EnclosureType> current_reach_enclosures;
        ListSet<EnclosureType> current_evolve_enclosures;

        make_lpair(current_reach_enclosures,current_evolve_enclosures) = evolver.reach_evolve(initial_enclosure,time,semantics);

        for(auto enclosure : current_reach_enclosures) {
            enclosure.adjoin_outer_approximation_to(reach_cells,accuracy);
        }
        ARIADNE_LOG_PRINTLN_AT(1,"final reach size = "<<reach_cells.size());

        for(auto enclosure : current_evolve_enclosures) {
            ARIADNE_LOG_PRINTLN_AT(2,"enclosure = "<<enclosure);
            enclosure.adjoin_outer_approximation_to(evolve_cells,accuracy);
            EnclosureType new_enclosure = enclosure;
            new_enclosure.clear_time();
            evolve_enclosures.adjoin(new_enclosure);
        }
    }
    ARIADNE_LOG_PRINTLN("final evolve size = "<<evolve_cells.size());
    ARIADNE_LOG_PRINTLN("final evolve enclosures size = "<<evolve_enclosures.size());
    return result;
}


// Helper functions for operators on lists of sets.
template<class SYS> auto
ReachabilityAnalyser<SYS>::
_upper_reach(const StorageType& set,
             const TimeType& time,
             const Nat accuracy,
             const EvolverType& evolver) const
    -> StorageType
{
    GridType grid=set.grid();
    StorageType result(grid,this->system());
    StorageType cells=set;
    cells.mince(accuracy);
    for(auto cell : cells) {
        EnclosureType initial_enclosure = evolver.enclosure(cell.box());
        ListSet<EnclosureType> reach = evolver.reach(initial_enclosure,time,Semantics::UPPER);
        for(auto enclosure : reach) {
            enclosure.adjoin_outer_approximation_to(result,accuracy);
        }
    }
    return result;
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
_upper_evolve(const StorageType& set,
              const TimeType& time,
              const Nat accuracy,
              const EvolverType& evolver) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    GridType grid=set.grid();
    StorageType result(grid,this->system()); StorageType cells=set; cells.mince(accuracy);
    for(auto cell : cells) {
        ARIADNE_LOG_PRINTLN_AT(1,"Evolving cell = "<<cell);
        EnclosureType initial_enclosure = evolver.enclosure(cell.box());
        ListSet<EnclosureType> final = evolver.evolve(initial_enclosure,time,Semantics::UPPER);
        for(auto enclosure : final) {
            enclosure.adjoin_outer_approximation_to(result,accuracy);
        }
    }
    ARIADNE_LOG_PRINTLN("_upper_evolve result size = "<<result.size());
    return result;
}


template<class SYS> Void
ReachabilityAnalyser<SYS>::
_adjoin_upper_reach_evolve(StorageType& reach_cells,
                           StorageType& evolve_cells,
                           const StorageType& set,
                           const TimeType& time,
                           const Nat accuracy,
                           const EvolverType& evolver) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    GridType grid=set.grid();
    StorageType cells=set;
    cells.mince(accuracy);

    ARIADNE_LOG_PRINTLN("Evolving "<<cells.size()<<" cells");
    ProgressIndicator indicator(cells.size());
    for(auto const& cell : cells) {
        ARIADNE_LOG_PRINTLN_AT(1,"evolving cell = "<<cell);
        EnclosureType initial_enclosure = evolver.enclosure(cell.box());
        ListSet<EnclosureType> reach_enclosures;
        ListSet<EnclosureType> final_enclosures;
        make_lpair(reach_enclosures,final_enclosures) = evolver.reach_evolve(initial_enclosure,time,Semantics::UPPER);
        ARIADNE_LOG_PRINTLN_AT(1,"computed "<<reach_enclosures.size()<<" reach enclosures and "<<final_enclosures.size()<<" final enclosures.");

        for(auto enclosure : reach_enclosures) {
            enclosure.adjoin_outer_approximation_to(reach_cells,accuracy);
        }
        for(auto enclosure : final_enclosures) {
            enclosure.adjoin_outer_approximation_to(evolve_cells,accuracy);
        }
        ARIADNE_LOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% of cells evolved.");
    }

    ARIADNE_LOG_PRINTLN("final reach size = "<<reach_cells.size());
    ARIADNE_LOG_PRINTLN("final evolve size = "<<evolve_cells.size());
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
lower_evolve(const OvertSetInterfaceType& initial_set,
             const TimeType& time) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    Nat grid_fineness = this->_configuration->maximum_grid_fineness();
    Nat grid_extent = this->_configuration->maximum_grid_extent();
    GridType grid=this->_configuration->grid();
    StorageType initial_cells(grid,this->system());
    StorageType final_cells(grid,this->system());

    ARIADNE_LOG_PRINTLN("Computing lower evolve set...");

    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_extent,grid_fineness+4);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_cells.size()="<<initial_cells.size());
    for(auto cell : initial_cells) {
        EnclosureType initial_enclosure=_evolver->enclosure(cell.box());
        ListSet<EnclosureType> final_enclosures=_evolver->evolve(initial_enclosure,time,Semantics::LOWER);
        for(auto enclosure : final_enclosures) {
            enclosure.adjoin_outer_approximation_to(final_cells,grid_fineness);
        }
    }
    return final_cells;
}



template<class SYS> auto
ReachabilityAnalyser<SYS>::
lower_reach(const OvertSetInterfaceType& initial_set,
            const TimeType& time) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    Nat grid_fineness = this->_configuration->maximum_grid_fineness();
    Nat grid_extent = this->_configuration->maximum_grid_extent();
    const GridType& grid=this->_configuration->grid();
    StorageType initial_cells(grid,this->system());
    StorageType reach_cells(grid,this->system());

    ARIADNE_LOG_PRINTLN("Computing lower reach set...");

    ARIADNE_LOG_PRINTLN_AT(1,"Adjoining initial set to the grid...");
    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_extent,grid_fineness+4);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_cells.size()="<<initial_cells.size());
    for(auto cell : initial_cells) {
        EnclosureType initial_enclosure=_evolver->enclosure(cell.box());
        ListSet<EnclosureType> reach_enclosures=_evolver->reach(initial_enclosure,time,Semantics::LOWER);
        for(auto enclosure : reach_enclosures) {
            enclosure.adjoin_outer_approximation_to(reach_cells,grid_fineness);
        }
    }
    return reach_cells;
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
lower_reach_evolve(const OvertSetInterfaceType& initial_set,
                   const TimeType& time) const
    -> Pair<StorageType,StorageType>
{
    ARIADNE_LOG_SCOPE_CREATE;
    Nat grid_fineness = this->_configuration->maximum_grid_fineness();
    Nat grid_extent = this->_configuration->maximum_grid_extent();

    ARIADNE_LOG_PRINTLN("Computing lower reach evolve...");

    const GridType& grid=this->_configuration->grid();

    StorageType initial_cells(grid,this->system());

    StorageType reach(grid,this->system()); StorageType evolve_cells(grid,this->system());

    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_extent,grid_fineness+4);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_cells.size()="<<initial_cells.size());

    for(auto cell : initial_cells) {
        EnclosureType initial_enclosure=_evolver->enclosure(cell.box());
        ListSet<EnclosureType> reach_enclosures;
        ListSet<EnclosureType> final_enclosures;
        make_lpair(reach_enclosures,final_enclosures) = _evolver->reach_evolve(initial_enclosure,time,Semantics::LOWER);
        for(auto enclosure : reach_enclosures) {
            enclosure.adjoin_outer_approximation_to(reach,grid_fineness);
        }
        for(auto enclosure : final_enclosures) {
            enclosure.adjoin_outer_approximation_to(evolve_cells,grid_fineness);
        }
    }
    return make_pair(reach,evolve_cells);
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
lower_reach(const OvertSetInterfaceType& initial_set) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    TimeType transient_time = this->_configuration->transient_time();
    TimeType lock_to_grid_time=this->_configuration->lock_to_grid_time();
    Nat maximum_grid_fineness = this->_configuration->maximum_grid_fineness();
    Nat maximum_grid_extent = this->_configuration->maximum_grid_extent();
    ARIADNE_LOG_PRINTLN_AT(1,"transient_time=("<<transient_time<<")");
    ARIADNE_LOG_PRINTLN_AT(1,"lock_to_grid_time=("<<lock_to_grid_time<<")");
    ARIADNE_LOG_PRINTLN_AT(1,"initial_set="<<initial_set);

    const GridType& grid=this->_configuration->grid();

    Bool has_bounding_domain = this->_configuration->bounding_domain_ptr() != nullptr;

    StorageType bounding(grid,this->system());
    if (has_bounding_domain)
        bounding.adjoin_outer_approximation(this->_configuration->bounding_domain(),maximum_grid_fineness);

    ARIADNE_LOG_PRINTLN_AT(1,"maximum_grid_extent="<<maximum_grid_extent);
    ARIADNE_LOG_PRINTLN_AT(1,"bounding_size="<<bounding.size());

    StorageType initial_cells(grid,this->system()), evolve_cells(grid,this->system());
    initial_cells.adjoin_lower_approximation(initial_set,maximum_grid_extent,maximum_grid_fineness+4);

    if (has_bounding_domain) initial_cells.restrict(bounding);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_size="<<initial_cells.size());

    ListSet<EnclosureType> starting_enclosures;

    StorageType reach_cells(grid,this->system());

    if(possibly(transient_time > TimeType(0))) {
        ARIADNE_LOG_PRINTLN("Computing transient evolution...");

        ListSet<EnclosureType> initial_enclosures;
        for(auto cell : initial_cells)
            initial_enclosures.adjoin(_evolver->enclosure(cell.box()));

        make_lpair(reach_cells,evolve_cells) =
            _reach_evolve_resume(initial_enclosures,transient_time,
                maximum_grid_fineness,starting_enclosures,Semantics::LOWER,*_evolver);

        evolve_cells.restrict_to_extent(maximum_grid_extent);
        ARIADNE_LOG_PRINTLN("Completed restriction to extent.");
        if (has_bounding_domain) reach_cells.restrict(bounding);

        ARIADNE_LOG_PRINTLN_AT(1,"reach size="<<reach_cells.size());
        ARIADNE_LOG_PRINTLN_AT(1,"evolve size="<<evolve_cells.size());
        ARIADNE_LOG_PRINTLN("Found "<<reach_cells.size()<<" cells.");
    }
    ARIADNE_LOG_PRINTLN("Computing recurrent evolution...");

    while(!evolve_cells.is_empty()) {

        StorageType new_reach_cells(grid,this->system());
        ListSet<EnclosureType> new_evolve_enclosures;
        make_lpair(new_reach_cells,evolve_cells) = _reach_evolve_resume(starting_enclosures,lock_to_grid_time,
                maximum_grid_fineness,new_evolve_enclosures,Semantics::LOWER,*_evolver);

        ARIADNE_LOG_PRINTLN_AT(1,"Removing already reached cells...");
        evolve_cells.remove(reach_cells);
        ARIADNE_LOG_PRINTLN_AT(1,"Restricting to extent...");
        evolve_cells.restrict_to_extent(maximum_grid_extent);
        ARIADNE_LOG_PRINTLN_AT(1,"Restricting to bounding...");
        if (has_bounding_domain) evolve_cells.restrict(bounding);
        ARIADNE_LOG_PRINTLN_AT(1,"Adjoining...");
        reach_cells.adjoin(new_reach_cells);
        ARIADNE_LOG_PRINTLN_AT(1,"Found "<<new_reach_cells.size()<<" cells, of which "<<evolve_cells.size()<<" are new.");

        starting_enclosures = new_evolve_enclosures;
    }
    reach_cells.recombine();
    reach_cells.restrict_to_extent(maximum_grid_extent);
    if (has_bounding_domain) reach_cells.restrict(bounding);

    return reach_cells;
}



inline Natural compute_time_steps(Real time, Real lock_to_grid_time) {
    return cast_positive(round(time/lock_to_grid_time));
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
upper_evolve(const CompactSetInterfaceType& initial_set,
             const TimeType& time) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    const GridType& grid=this->_configuration->grid();
    StorageType evolve_cells(grid,this->system());
    Nat grid_fineness = this->_configuration->maximum_grid_fineness();
    evolve_cells.adjoin_outer_approximation(initial_set,grid_fineness);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_evolve.size()="<<evolve_cells.size());
    TimeType lock_to_grid_time=this->_configuration->lock_to_grid_time();
    Natural time_steps=compute_time_steps(time,lock_to_grid_time);
    TimeType remainder_time=time-(time_steps*lock_to_grid_time);
    ARIADNE_LOG_PRINTLN_AT(1,"real_time="<<time);
    ARIADNE_LOG_PRINTLN_AT(1,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time);

    for(Natural i=0u; i!=time_steps; ++i) {
        ARIADNE_LOG_PRINTLN_AT(2,"computing "<<i+1u<<"-th reachability step...");
        evolve_cells=this->_upper_evolve(evolve_cells,lock_to_grid_time,grid_fineness,*_evolver);
    }
    ARIADNE_LOG_PRINTLN_AT(1,"remainder_time="<<remainder_time);
    if(!evolve_cells.is_empty() && possibly(remainder_time > 0)) {
        ARIADNE_LOG_PRINTLN_AT(1,"computing evolution for remainder time...");
        evolve_cells=this->_upper_evolve(evolve_cells,remainder_time,grid_fineness,*_evolver);
    }
    evolve_cells.recombine();
    ARIADNE_LOG_PRINTLN("final_evolve.size()="<<evolve_cells.size());
    return evolve_cells;
}



template<class SYS> auto
ReachabilityAnalyser<SYS>::
upper_reach(const CompactSetInterfaceType& initial_set,
            const TimeType& time) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"initial_set="<<initial_set);
    const GridType& grid=this->_configuration->grid();
    StorageType evolve_cells(grid,this->system());
    Nat grid_fineness = this->_configuration->maximum_grid_fineness();
    ARIADNE_LOG_PRINTLN_AT(2,"grid_fineness="<<grid_fineness);
    evolve_cells.adjoin_outer_approximation(initial_set,grid_fineness);
    ARIADNE_LOG_PRINTLN("initial size = "<<evolve_cells.size());
    StorageType reach_cells(evolve_cells);
    ARIADNE_LOG_PRINTLN("reach size ="<<reach_cells.size());
    TimeType lock_to_grid_time = this->_configuration->lock_to_grid_time();
    Natural time_steps=compute_time_steps(time,lock_to_grid_time);
    TimeType remainder_time=time-time_steps*lock_to_grid_time;

    ARIADNE_LOG_PRINTLN_AT(1,"time="<<time);
    ARIADNE_LOG_PRINTLN_AT(1,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time);
    StorageType found_cells(grid,this->system());
    StorageType accumulated_evolve_cells(grid,this->system());
    ARIADNE_LOG_PRINTLN("Computing time steps...")

    for(Natural i=0u; i!=time_steps; ++i) {
        accumulated_evolve_cells.adjoin(found_cells);
        ARIADNE_LOG_PRINTLN_AT(1,"computing "<<i+1u<<"-th reachability step...");
        this->_adjoin_upper_reach_evolve(found_cells,evolve_cells,evolve_cells,lock_to_grid_time,grid_fineness,*_evolver);
        ARIADNE_LOG_PRINTLN_AT(1,"found.size()="<<found_cells.size());
        ARIADNE_LOG_PRINTLN_AT(1,"evolve.size()="<<evolve_cells.size());
        evolve_cells.remove(accumulated_evolve_cells);
        reach_cells.adjoin(found_cells);
        accumulated_evolve_cells.adjoin(evolve_cells);
        ARIADNE_LOG_PRINTLN_AT(1,"found "<<found_cells.size()<<" cells, with "<<evolve_cells.size()<<" new intermediate.");
        if(evolve_cells.is_empty()) break;
        ARIADNE_LOG_PRINTLN_AT(1,"evolve_cells="<<evolve_cells);
    }
    ARIADNE_LOG_PRINTLN("remainder_time="<<remainder_time);
    if(!evolve_cells.is_empty() && possibly(remainder_time > 0)) {
        ARIADNE_LOG_PRINTLN("computing evolution for remainder time...");
        this->_adjoin_upper_reach_evolve(found_cells,accumulated_evolve_cells,evolve_cells,remainder_time,grid_fineness,*_evolver);
        reach_cells.adjoin(found_cells);
    }
    // This last step is necessary to add the final set to the result.
    reach_cells.adjoin(evolve_cells);
    reach_cells.recombine();
    ARIADNE_LOG_PRINTLN("final_reach size = "<<reach_cells.size());
    return reach_cells;
}



template<class SYS> auto
ReachabilityAnalyser<SYS>::
upper_reach_evolve(const CompactSetInterfaceType& initial_set,
                   const TimeType& time) const
    -> Pair<StorageType,StorageType>
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"initial_set="<<initial_set);
    const GridType& grid=this->_configuration->grid();
    StorageType evolve_cells(grid,this->system());
    Nat grid_fineness = this->_configuration->maximum_grid_fineness();
    ARIADNE_LOG_PRINTLN_AT(1,"grid_fineness="<<grid_fineness);
    evolve_cells.adjoin_outer_approximation(initial_set,grid_fineness);

    ARIADNE_LOG_PRINTLN_AT(1,"initial_evolve"<<evolve_cells);
    StorageType reach_cells(evolve_cells);
    ARIADNE_LOG_PRINTLN_AT(1,"reach="<<reach_cells);
    TimeType lock_to_grid_time = this->_configuration->lock_to_grid_time();
    Natural time_steps=compute_time_steps(time,lock_to_grid_time);
    TimeType remainder_time=time-time_steps*lock_to_grid_time;
    ARIADNE_LOG_PRINTLN("time="<<time);
    ARIADNE_LOG_PRINTLN("time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time);

    ARIADNE_LOG_PRINTLN("Computing reachability steps...");

    StorageType found_cells(grid,this->system());
    for(Natural i=0u; i!=time_steps; ++i) {
        ARIADNE_LOG_PRINTLN_AT(1,"computing "<<i+1u<<"-th reachability step...");
        this->_adjoin_upper_reach_evolve(found_cells,evolve_cells,evolve_cells,lock_to_grid_time,grid_fineness,*_evolver);
        ARIADNE_LOG_PRINTLN_AT(2,"evolve.size()="<<evolve_cells.size());
        reach_cells.adjoin(found_cells);
        ARIADNE_LOG_PRINTLN_AT(1,"found "<<found_cells.size()<<" cells.");
    }
    ARIADNE_LOG_PRINTLN("remainder_time="<<remainder_time);
    if(!evolve_cells.is_empty() && possibly(remainder_time > 0)) {
        ARIADNE_LOG_PRINTLN("computing evolution for remainder time...");
        this->_adjoin_upper_reach_evolve(found_cells,evolve_cells,evolve_cells,remainder_time,grid_fineness,*_evolver);
        reach_cells.adjoin(found_cells);
    }
    reach_cells.recombine();
    ARIADNE_LOG_PRINTLN("reach="<<reach_cells);
    evolve_cells.recombine();
    ARIADNE_LOG_PRINTLN("evolve="<<evolve_cells);
    return std::make_pair(reach_cells,evolve_cells);
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
outer_chain_reach(const CompactSetInterfaceType& initial_set) const
    -> StorageType
{
    ARIADNE_LOG_SCOPE_CREATE;
    TimeType transient_time = this->_configuration->transient_time();
    TimeType lock_to_grid_time=this->_configuration->lock_to_grid_time();
    Nat maximum_grid_fineness = this->_configuration->maximum_grid_fineness();
    Nat maximum_grid_extent = this->_configuration->maximum_grid_extent();
    ARIADNE_LOG_PRINTLN_AT(1,"transient_time=("<<transient_time<<")");
    ARIADNE_LOG_PRINTLN_AT(1,"lock_to_grid_time=("<<lock_to_grid_time<<")");
    ARIADNE_LOG_PRINTLN_AT(1,"initial_set="<<initial_set);

    const GridType& grid=this->_configuration->grid();

    Bool has_bounding_domain = this->_configuration->bounding_domain_ptr()!=nullptr;

    StorageType bounding(grid,this->system());
    if (has_bounding_domain)
        bounding.adjoin_outer_approximation(this->_configuration->bounding_domain(),maximum_grid_fineness);

    ARIADNE_LOG_PRINTLN_AT(1,"maximum_grid_extent="<<maximum_grid_extent);
    ARIADNE_LOG_PRINTLN_AT(1,"bounding_size="<<bounding.size());

    StorageType initial_cells(grid,this->system());
    initial_cells.adjoin_outer_approximation(initial_set,maximum_grid_fineness);
    _checked_restriction(initial_cells,bounding);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_size="<<initial_cells.size());

    Nat stage=0;

    StorageType reach_cells(grid,this->system());
    StorageType evolve_cells(grid,this->system());
    if(definitely(transient_time > TimeType(0))) {
        ARIADNE_LOG_PRINTLN("Computing transient evolution...");
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,initial_cells,transient_time,maximum_grid_fineness,*_evolver);
        _checked_restriction(evolve_cells,bounding);
        evolve_cells.recombine();
        evolve_cells.mince(maximum_grid_fineness);
        ARIADNE_LOG_PRINTLN_AT(1,"transient evolve_size="<<evolve_cells.size());
        ARIADNE_LOG_PRINTLN("Found "<<reach_cells.size()<<" cells.");
    } else {
        evolve_cells=initial_cells;
    }

    ARIADNE_LOG_PRINTLN("Computing recurrent evolution...");
    StorageType starting_cells = evolve_cells;
    StorageType accumulated_evolve_cells = evolve_cells;

    while(!starting_cells.is_empty()) {
        ++stage;
        ARIADNE_LOG_PRINTLN_AT(1,"stage="<<std::setw(3)<<stage<<
                      " #starting="<<std::setw(4)<<std::left<<starting_cells.size()<<
                      " #reached="<<std::setw(4)<<std::left<<reach_cells.size()<<
                      " #evolved="<<std::setw(4)<<std::left<<accumulated_evolve_cells.size());
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,starting_cells,
                                         lock_to_grid_time,maximum_grid_fineness,*_evolver);
        ARIADNE_LOG_PRINTLN_AT(2,"reach.size()="<<reach_cells.size());
        ARIADNE_LOG_PRINTLN_AT(2,"evolve.size()="<<evolve_cells.size());
        _checked_restriction(evolve_cells,bounding);
        starting_cells = evolve_cells;
        starting_cells.remove(accumulated_evolve_cells);
        accumulated_evolve_cells.adjoin(starting_cells);
        starting_cells.mince(maximum_grid_fineness);
        ARIADNE_LOG_PRINTLN_AT(2,"evolved to "<<evolve_cells.size()<<" cells, of which "<<starting_cells.size()<<" are new.");
    }
    _checked_restriction(reach_cells,bounding);
    return reach_cells;
}


template<class SYS> auto
ReachabilityAnalyser<SYS>::
verify_safety(const CompactSetInterfaceType& initial_set,
              const OpenSetInterfaceType& safe_set) const
    -> SafetyCertificateType
{
//    const BoundedSetInterfaceType* bounded_safe_set_ptr=dynamic_cast<BoundedSetInterfaceType const*>(&safe_set);
//    assert(bounded_safe_set_ptr != nullptr);
//
//    BoundingDomainType safe_set_bounding_box=cast_exact(bounded_safe_set_ptr->bounding_box());
//    const GridType& grid=this->_base_configuration->grid();
//    StorageType safe_cells=inner_approximation(safe_set, grid, safe_set_bounding_box, this->_base_configuration->maximum_grid_fineness());

    ARIADNE_LOG_SCOPE_CREATE;
    const RegularLocatedSetInterfaceType* bounded_safe_set_ptr=dynamic_cast<RegularLocatedSetInterfaceType const*>(&safe_set);
    assert(bounded_safe_set_ptr != nullptr);
    const GridType& grid=this->_configuration->grid();
    auto safe_paving=inner_approximation(*bounded_safe_set_ptr, grid, this->_configuration->maximum_grid_fineness());
    StorageType safe_cells(safe_paving,this->system());

    TimeType transient_time = this->_configuration->transient_time();
    TimeType lock_to_grid_time=this->_configuration->lock_to_grid_time();
    Nat maximum_grid_fineness = this->_configuration->maximum_grid_fineness();
    ARIADNE_LOG_PRINTLN_AT(1,"transient_time=("<<transient_time<<")");
    ARIADNE_LOG_PRINTLN_AT(1,"lock_to_grid_time=("<<lock_to_grid_time<<")");
    ARIADNE_LOG_PRINTLN_AT(1,"initial_set="<<initial_set);

    StorageType initial_cells(grid,this->system());
    initial_cells.adjoin_outer_approximation(initial_set,maximum_grid_fineness);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_size="<<initial_cells.size());

    if(not subset(initial_cells,safe_cells)) {
        return SafetyCertificateType { indeterminate,initial_cells,safe_cells };
    }

    StorageType reach_cells(grid,this->system());
    StorageType evolve_cells(grid,this->system());

    if(definitely(transient_time > 0)) {
        ARIADNE_LOG_PRINTLN("Computing transient evolution...");
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,initial_cells,transient_time,maximum_grid_fineness,*_evolver);
        evolve_cells.mince(maximum_grid_fineness);
        ARIADNE_LOG_PRINTLN_AT(1,"transient_reach_size="<<reach_cells.size());
        ARIADNE_LOG_PRINTLN_AT(1,"transient evolve_size="<<evolve_cells.size());
        ARIADNE_LOG_PRINTLN_AT(1,"found "<<reach_cells.size()<<" cells.");
    } else {
        reach_cells = initial_cells;
        evolve_cells = initial_cells;
    }

    ARIADNE_LOG_PRINTLN("Computing recurrent evolution...");
    StorageType starting_cells = evolve_cells;
    StorageType current_evolve_cells(grid,this->system());

    Nat stage=0;
    while(!starting_cells.is_empty()) {
        ARIADNE_LOG_PRINTLN_AT(1,"i="<<std::setw(3)<<++stage<<" #s="<<std::setw(4)<<std::left<<starting_cells.size());
        if(not subset(reach_cells,safe_cells)) { return SafetyCertificateType { indeterminate,reach_cells,safe_cells }; }

        current_evolve_cells = evolve_cells;
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,starting_cells,
                                         lock_to_grid_time,maximum_grid_fineness,*_evolver);
        evolve_cells.mince(maximum_grid_fineness);
        ARIADNE_LOG_PRINTLN_AT(2,"reach.size()="<<reach_cells.size());
        ARIADNE_LOG_PRINTLN_AT(2,"evolve.size()="<<evolve_cells.size());
        starting_cells = evolve_cells;
        starting_cells.remove(current_evolve_cells);
        starting_cells.mince(maximum_grid_fineness);
        ARIADNE_LOG_PRINTLN_AT(2,"evolved to "<<evolve_cells.size()<<" cells, of which "<<starting_cells.size()<<" are new.");
    }

    return SafetyCertificateType { true,reach_cells,safe_cells };
}


template<class SYS> Void
ReachabilityAnalyser<SYS>::
_checked_restriction(StorageType& set, const StorageType& bounding) const
{
    ChainOverspillPolicy policy = this->_configuration->outer_overspill_policy();
    StorageType set_copy(set.grid(),set.auxiliary_data());

    if (policy != ChainOverspillPolicy::IGNORE) {
        set_copy = set;
    }
    set.restrict_to_extent(this->_configuration->maximum_grid_extent());
    if (!bounding.is_empty()) set.restrict(bounding);

    if (policy != ChainOverspillPolicy::IGNORE) {
        set_copy.remove(set);
        if (!set_copy.is_empty()) {
            if (policy == ChainOverspillPolicy::WARNING) {
                ARIADNE_WARN("The computed chain reach has been restricted, an outer approximation is .");
            }
            if (policy == ChainOverspillPolicy::ERROR) {
                ARIADNE_THROW(OuterChainOverspill,"outer_chain_reach","The computed chain reach has been restricted.");
            }
        }
    }
}

template<class SYS> Void ReachabilityAnalyserConfiguration<SYS>::set_grid(const GridType& grid) {
    this->_grid_ptr.reset(new GridType(grid));
}

template<class SYS> Void ReachabilityAnalyserConfiguration<SYS>::set_bounding_domain(const BoundingDomainType& bounding_domain) {
    this->_bounding_domain_ptr.reset(new BoundingDomainType(bounding_domain));
}


template<class SYS> ReachabilityAnalyserConfiguration<SYS>::ReachabilityAnalyserConfiguration(ReachabilityAnalyser<SYS>& analyser)
    : _analyser(analyser)
{
    set_transient_time(0);
    set_lock_to_grid_time(1);
    set_maximum_grid_fineness(3);
    set_maximum_grid_extent(16);
    set_bounding_domain(BoundingDomainType(_analyser.system().dimension(),{-1,+1}));
    set_grid(GridType(_analyser.system()));
    set_outer_overspill_policy(ChainOverspillPolicy::ERROR);
}


template<class SYS>
OutputStream&
ReachabilityAnalyserConfiguration<SYS>::_write(OutputStream& os) const
{
    os << "ReachabilityAnalyserSettings"
       << "(\n  transient_time=" << transient_time()
       << ",\n  lock_to_grid_time=" << lock_to_grid_time()
       << ",\n  maximum_grid_fineness=" << maximum_grid_fineness()
       << ",\n  maximum_grid_extent=" << maximum_grid_extent()
       << ",\n  bounding_domain=" << bounding_domain()
       << ",\n  grid=" << grid()
       << "\n)\n";
    return os;
}


} // namespace Ariadne
