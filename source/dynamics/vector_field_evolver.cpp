/***************************************************************************
 *            dynamics/vector_field_evolver.cpp
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
#include "utility/container.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"

#include "solvers/integrator.hpp"

#include "output/logging.hpp"
#include "output/progress_indicator.hpp"

#include "dynamics/vector_field.hpp"
#include "dynamics/vector_field_evolver.hpp"

#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"

namespace Ariadne {

#warning
inline PositiveValidatedUpperNumber abs(PositiveValidatedUpperNumber y) { return y; }

namespace {



template<class ES> List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES,ES> split=enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

} // namespace

EffectiveVectorMultivariateFunction make_auxiliary_function(
    Space<Real> const& state_space,
    List<RealAssignment> const& algebraic);

EffectiveVectorMultivariateFunction make_dynamic_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<DottedRealAssignment> const& differential);


VectorField::VectorField(List<DottedRealAssignment> const& dynamics)
    : VectorField(dynamics, List<RealAssignment>())
{
}

VectorField::VectorField(List<DottedRealAssignment> const& dynamics, List<RealAssignment> const& auxiliary)
    : _dynamics(dynamics), _auxiliary(auxiliary)
    , _dynamic_function(make_dynamic_function(left_hand_sides(dynamics),auxiliary,dynamics))
    , _auxiliary_function(make_auxiliary_function(left_hand_sides(dynamics),auxiliary))
{
}

VectorField::VectorField(EffectiveVectorMultivariateFunction const& function)
    : _dynamic_function(function), _auxiliary_function(0u,function.domain())
{
    ARIADNE_PRECONDITION(function.result_size()==function.argument_size());
}

RealSpace VectorField::state_space() const {
    return RealSpace(left_hand_sides(this->_dynamics));
}

RealSpace VectorField::auxiliary_space() const {
    return RealSpace(left_hand_sides(this->_auxiliary));
}


OutputStream& operator<<(OutputStream& os, const VectorField& vf) {
    os << "VectorField( dynamic_function = " << vf.dynamic_function() << ", "
          "auxiliary_function = " << vf.auxiliary_function() << ", "
          "dynamics = " << vf._dynamics << ", "
          "auxiliary = " << vf._auxiliary << ")";
    return os;
}

// Allow subdivisions in upper evolution
const Bool ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution
const Bool ENABLE_PREMATURE_TERMINATION = false;

using std::shared_ptr;

class DegenerateCrossingException { };



VectorFieldEvolver::VectorFieldEvolver(const SystemType& system, const IntegratorInterface& i)
    : _sys_ptr(system.clone())
    , _integrator(i.clone())
    , _configuration(new ConfigurationType())
{
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const ExactBoxType& box) const {
    return EnclosureType(box,this->system().state_space(),this->function_factory());
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealBox& box) const {
    return EnclosureType(box,this->system().state_space(),this->function_factory());
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealVariablesBox& box) const {
    return EnclosureType(box,this->system().state_space(),this->function_factory());
}

#warning
//typename VectorFieldEvolver::FunctionFactoryType const& VectorFieldEvolver::function_factory() const {
//    return std::dynamic_pointer_cast<const IntegratorBase>(this->_integrator)->function_factory();
//}


Orbit<VectorFieldEvolver::EnclosureType>
VectorFieldEvolver::
orbit(const EnclosureType& initial_set,
      const TimeType& time,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     initial_set,time,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}

Void VectorFieldEvolver::
_append_initial_set(List<TimedEnclosureType>& working_sets,
                   const TimeStepType& initial_time,
                   const EnclosureType& current_set) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    if (possibly(current_set.euclidean_set().bounding_box().radius() > this->_configuration->maximum_enclosure_radius())) {
        ARIADNE_LOG_PRINTLN("initial set too large, splitting");
        Pair<EnclosureType,EnclosureType> split_sets = current_set.split();
        if(!definitely(split_sets.first.is_empty())) { _append_initial_set(working_sets,initial_time,split_sets.first); }
        if(!definitely(split_sets.second.is_empty())) { _append_initial_set(working_sets,initial_time,split_sets.second); }
    } else {
        working_sets.push_back(make_pair(initial_time,current_set));
    }
}

Void test_print() {
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("test\r\n");
}


Void
VectorFieldEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const EnclosureType& initial_set,
           const TimeType& maximum_time,
           Semantics semantics,
           Bool reach) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    List< TimedEnclosureType > working_sets;

    {
        TimeStepType initial_time = 0u;
        // Append the initial set, possibly splitting it
        _append_initial_set(working_sets,initial_time,initial_set);
    }

    ProgressIndicator initials_indicator(working_sets.size());
    ProgressIndicator time_indicator(maximum_time.get_d());

    while(!working_sets.empty()) {
        TimedEnclosureType current_timed_set=working_sets.back();
        working_sets.pop_back();
        TimeStepType current_time=current_timed_set.first;
        EnclosureType current_set_model=current_timed_set.second;
        FloatDPUpperBound current_set_radius=current_set_model.euclidean_set().bounding_box().radius();
        if(definitely(current_time>=maximum_time)) {
            final_sets.adjoin(current_set_model);
        } else if(semantics == Semantics::UPPER && ENABLE_SUBDIVISIONS
                  && decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List< EnclosureType > subdivisions=subdivide(current_set_model);
            for(SizeType i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_pair(current_time,subdivided_set_model));
            }
        } else if(semantics == Semantics::LOWER && ENABLE_PREMATURE_TERMINATION && decide(current_set_radius>this->_configuration->maximum_enclosure_radius())) {
            ARIADNE_WARN("Terminating lower evolution at time " << current_time << " and set " << current_set_model << " due to maximum radius being exceeded.")
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  current_timed_set,maximum_time,
                                  semantics,reach);
        }

        ARIADNE_LOG_PRINTLN("#w="<<std::setw(4)<<working_sets.size()
                            <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                            <<" t="<<std::setw(7)<<std::fixed<<current_time.get_d()
                            <<" p="<<std::setw(4)<<std::left<<current_set_model.number_of_parameters()
                            <<" r="<<std::setw(7)<<current_set_model.radius()
                            <<" c="<<current_set_model.centre());

        initials_indicator.update_current(final_sets.size());
        time_indicator.update_current(current_time.get_d());
        if (initials_indicator.final_value() > 1) { ARIADNE_LOG_SCOPE_PRINTHOLD("[" << time_indicator.symbol() << "] " << initials_indicator.percentage() << "% of sets, " << time_indicator.percentage() << "% of current set"); }
        else ARIADNE_LOG_SCOPE_PRINTHOLD("[" << time_indicator.symbol() << "] " << time_indicator.percentage() << "%");
    }
}





Void
VectorFieldEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& working_timed_set_model,
                const TimeType& maximum_time,
                Semantics semantics,
                Bool reach) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    typedef EffectiveVectorMultivariateFunction FunctionType;

    EnclosureType current_set_model;
    TimeStepType current_time;
    ARIADNE_LOG_PRINTLN_AT(1,"working_timed_set_model = "<<working_timed_set_model);
    make_lpair(current_time,current_set_model)=working_timed_set_model;

    ARIADNE_LOG_PRINTLN("current_time = "<<current_time);
    ARIADNE_LOG_PRINTLN("current_set_model = "<<current_set_model);

    ARIADNE_LOG_PRINTLN("box = "<<current_set_model.bounding_box());
    ARIADNE_LOG_PRINTLN("radius = "<<current_set_model.euclidean_set().bounding_box().radius());
    //const SizeType nd=initial_set_model.result_size();
    //const SizeType ng=initial_set_model.argument_size();


    // Test to see if set requires reconditioning
    if(this->_configuration->enable_reconditioning() &&
       possibly(norm(current_set_model.state_function().errors()) > this->_configuration->maximum_spacial_error())) {

        ARIADNE_LOG_PRINTLN("reconditioning from errors "<<current_set_model.state_function().errors());
        current_set_model.recondition();
        working_sets.append(make_pair(current_time,current_set_model));
        return;
    }

    /////////////// Main Evolution ////////////////////////////////
    const FunctionType& dynamic=_sys_ptr->function();

    // Set evolution parameters
    const StepSizeType maximum_step_size=this->_configuration->maximum_step_size();
    //const FloatDP maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
    //const FloatDP zero_time=0.0;

    // Get bounding boxes for time and space bounding_box
    auto current_set_bounds=cast_exact_box(current_set_model.euclidean_set().bounding_box());
    ARIADNE_LOG_PRINTLN("current_set_bounds = "<<current_set_bounds);


    // Compute flow model
    // TODO: Modify this for general integrator interface
    //TaylorPicardIntegrator const* taylor_integrator=dynamic_cast<const TaylorPicardIntegrator*>(this->_integrator.operator->());
    IntegratorInterface const* integrator=this->_integrator.operator->();
    StepSizeType step_size=maximum_step_size;
    FlowStepModelType flow_model=integrator->flow_step(dynamic,current_set_bounds,step_size);
    ARIADNE_LOG_PRINTLN("step_size = "<<step_size);
    ARIADNE_LOG_PRINTLN_AT(1,"flow_model = "<<flow_model);
    FlowStepModelType flow_step_model=partial_evaluate(flow_model,flow_model.domain().size()-1u,step_size);
    ARIADNE_LOG_PRINTLN_AT(1,"flow_step_model = "<<flow_step_model);

    // Compute the integration time model
    TimeStepType next_time=current_time+TimeStepType(step_size);
    ARIADNE_LOG_PRINTLN_AT(1,"next_time = "<<next_time);
    // Compute the flow tube (reachable set) model and the final set
    EnclosureType reach_set_model=current_set_model;
    reach_set_model.apply_full_reach_step(flow_model);
    ARIADNE_LOG_PRINTLN_AT(1,"reach_set_model = "<<reach_set_model);
    EnclosureType next_set_model=current_set_model;
    next_set_model.apply_fixed_evolve_step(flow_model,step_size);
    ARIADNE_LOG_PRINTLN_AT(1,"next_set_model = "<<next_set_model);

    reach_sets.adjoin(reach_set_model);

    intermediate_sets.adjoin(EnclosureType(next_set_model));
    working_sets.push_back(make_pair(next_time,next_set_model));
}


VectorFieldEvolverConfiguration::VectorFieldEvolverConfiguration()
{
    set_maximum_step_size(1);
    set_maximum_enclosure_radius(100.0);
    set_enable_reconditioning(true);
    set_maximum_spacial_error(1e-2);
}


OutputStream&
VectorFieldEvolverConfiguration::_write(OutputStream& os) const
{
    os << "VectorFieldEvolverConfiguration("
       << "\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  enable_reconditioning=" << enable_reconditioning()
       << ",\n  maximum_spacial_error=" << maximum_spacial_error()
       << "\n)";
    return os;
}


}  // namespace Ariadne

