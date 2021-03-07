/***************************************************************************
 *            dynamics/map_evolver.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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
#include "config.hpp"

#include "utility/macros.hpp"
#include "utility/array.hpp"
#include "utility/tuple.hpp"
#include "utility/stlio.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "symbolic/assignment.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"

#include "output/logging.hpp"

#include "dynamics/map.hpp"
#include "dynamics/map_evolver.hpp"

namespace Ariadne {

namespace {

template<class ES> List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES,ES> split=enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

} // namespace


// The maximum allowable radius of a basic set.
MapEvolverConfiguration::RealType DEFAULT_MAXIMUM_ENCLOSURE_RADIUS = 100;
// Allow subdivisions in upper evolution.
const Bool DEFAULT_ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution.
const Bool DEFAULT_ENABLE_PREMATURE_TERMINATION = false;

using std::shared_ptr;

FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory();

class DegenerateCrossingException { };


EffectiveVectorMultivariateFunction make_auxiliary_function(
    Space<Real> const& state_space,
    List<RealAssignment> const& algebraic);

EffectiveVectorMultivariateFunction make_reset_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<PrimedRealAssignment> const& differential);


IteratedMap::IteratedMap(const EffectiveVectorMultivariateFunction& f)
    : _update_function(f)
{
    ARIADNE_PRECONDITION(f.result_size()==f.argument_size());
}

IteratedMap::IteratedMap(const List<PrimedRealAssignment>& updates)
    : IteratedMap(updates,List<RealAssignment>()) { }

IteratedMap::IteratedMap(const List<PrimedRealAssignment>& updates, List<RealAssignment> const& auxiliary)
    : _updates(updates), _auxiliary(auxiliary)
    , _update_function(make_reset_function(left_hand_sides(updates),auxiliary,updates))
    , _auxiliary_function(make_auxiliary_function(left_hand_sides(updates),auxiliary))
{
}

RealSpace IteratedMap::state_space() const {
    return RealSpace(left_hand_sides(this->_updates));
}

RealSpace IteratedMap::auxiliary_space() const {
    return RealSpace(left_hand_sides(this->_auxiliary));
}


OutputStream& operator<<(OutputStream& os, const IteratedMap& map) {
    os << "IteratedMap( update_function = " << map.update_function() << ", "
          "auxiliary_function = " << map.auxiliary_function() << ", "
          "updates = " << map._updates << ", "
          "auxiliary = " << map._auxiliary << ")";
    return os;
}




MapEvolver::MapEvolver(const SystemType& system)
    : _sys_ptr(system.clone())
    , _configuration(new ConfigurationType())
{
}

typename MapEvolver::FunctionFactoryType const MapEvolver::function_factory() const {
    return ValidatedFunctionModelDPFactory(make_taylor_function_factory());
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const ExactBoxType& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const ExactBoxType& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(), config);
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const RealBox& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const RealBox& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(), config);
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const RealVariablesBox& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(this->function_factory()));
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const RealVariablesBox& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(), config);
}

Orbit<MapEvolver::EnclosureType>
MapEvolver::
orbit(const EnclosureType& initial_set,
      const TerminationType& termination,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     initial_set,termination,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}

Void
MapEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const EnclosureType& initial_set,
           const TerminationType& maximum_time,
           Semantics semantics,
           Bool reach) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    List< TimedEnclosureType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG_PRINTLN_AT(1,"initial_set = "<<initial_set);
        EnclosureType initial_enclosure(initial_set);
        ARIADNE_LOG_PRINTLN_AT(1,"initial_enclosure = "<<initial_enclosure);
        TimeType initial_time = 0;
        ARIADNE_LOG_PRINTLN_AT(1,"initial_time = "<<initial_time);
        working_sets.push_back(make_pair(initial_time,initial_enclosure));
    }


    while(!working_sets.empty()) {
        TimedEnclosureType current_set=working_sets.back();
        working_sets.pop_back();
        EnclosureType initial_enclosure=current_set.second;
        TimeType initial_time=current_set.first;
        FloatDPUpperBound initial_set_radius=initial_enclosure.euclidean_set().bounding_box().radius();
        if(initial_time>=maximum_time) {
            final_sets.adjoin(EnclosureType(initial_enclosure));
        } else if(semantics == Semantics::UPPER && this->_configuration->enable_subdivisions()
                  && decide(initial_set_radius>this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List<EnclosureType> subdivisions=subdivide(initial_enclosure);
            for(SizeType i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_enclosure=subdivisions[i];
                working_sets.push_back(make_pair(initial_time,subdivided_enclosure));
            }
        } else if(semantics == Semantics::LOWER && this->_configuration->enable_premature_termination() && decide(initial_set_radius>this->_configuration->maximum_enclosure_radius())) {
            ARIADNE_WARN("Terminating lower evolution at time " << initial_time
                         << " and set " << initial_enclosure << " due to maximum radius being exceeded.");
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  current_set,maximum_time,
                                  semantics,reach);
        }
    }

}





Void
MapEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& current_set,
                const TimeType& maximum_time,
                Semantics semantics,
                Bool reach) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    EnclosureType initial_enclosure;
    TimeType initial_time;

    ARIADNE_LOG_PRINTLN_AT(1,"working_set = "<<current_set);
    make_lpair(initial_time,initial_enclosure)=current_set;

    ARIADNE_LOG_PRINTLN_AT(1,"initial_time = "<<initial_time);
    ARIADNE_LOG_PRINTLN_AT(1,"initial_enclosure = "<<initial_enclosure);

    ARIADNE_LOG_PRINTLN_AT(1,"box = "<<initial_enclosure.bounding_box()<<" ");
    ARIADNE_LOG_PRINTLN_AT(1,"radius = "<<initial_enclosure.euclidean_set().bounding_box().radius()<<"\n");
    //const SizeType nd=initial_enclosure.result_size();
    //const SizeType ng=initial_enclosure.argument_size();


    /////////////// Main Evolution ////////////////////////////////


    // Compute the map model
    EnclosureType& final_enclosure=initial_enclosure;
    final_enclosure.apply_map(_sys_ptr->function());
    TimeType final_time=initial_time+1;
    ARIADNE_LOG_PRINTLN("final_enclosure = "<<final_enclosure);

    reach_sets.adjoin(final_enclosure);

    intermediate_sets.adjoin(final_enclosure);
    working_sets.push_back(make_pair(final_time,final_enclosure));

}


MapEvolverConfiguration::MapEvolverConfiguration()
{
    set_maximum_enclosure_radius(DEFAULT_MAXIMUM_ENCLOSURE_RADIUS);
    set_enable_subdivisions(DEFAULT_ENABLE_SUBDIVISIONS);
    set_enable_premature_termination(DEFAULT_ENABLE_PREMATURE_TERMINATION);
}


OutputStream&
MapEvolverConfiguration::_write(OutputStream& os) const
{
    os << "MapEvolverSettings"
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << "\n)\n";
    return os;
}


}  // namespace Ariadne

