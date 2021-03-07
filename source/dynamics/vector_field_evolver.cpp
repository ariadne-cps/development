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
#include "function/taylor_function.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"

#include "solvers/integrator.hpp"

#include "output/logging.hpp"
#include "output/progress_indicator.hpp"

#include "concurrency/concurrency_manager.hpp"
#include "configuration/configuration_property.tpl.hpp"

#include "dynamics/vector_field.hpp"
#include "dynamics/vector_field_evolver.hpp"

#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"

namespace Ariadne {

template class RangeConfigurationProperty<ExactDouble>;
template class InterfaceListConfigurationProperty<IntegratorInterface>;

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

VectorFieldEvolver::VectorFieldEvolver(const SystemType& system, const ConfigurationType& configuration) :
    TaskRunnable(configuration), _sys_ptr(system.clone())
{ }

VectorFieldEvolver* VectorFieldEvolver::clone() const {
    return new VectorFieldEvolver(system(),configuration());
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const ExactBoxType& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(TaylorFunctionFactory(Sweeper<FloatDP>())));
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const ExactBoxType& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(),config);
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealBox& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(TaylorFunctionFactory(Sweeper<FloatDP>())));
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealBox& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(),config);
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealVariablesBox& box) const {
    return EnclosureType(box,this->system().state_space(),EnclosureConfiguration(TaylorFunctionFactory(Sweeper<FloatDP>())));
}

typename VectorFieldEvolver::EnclosureType VectorFieldEvolver::enclosure(const RealVariablesBox& box, const EnclosureConfiguration& config) const {
    return EnclosureType(box,this->system().state_space(),config);
}


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
    if (possibly(current_set.euclidean_set().bounding_box().radius() > configuration().maximum_enclosure_radius())) {
        ARIADNE_LOG_PRINTLN("initial set too large, splitting");
        Pair<EnclosureType,EnclosureType> split_sets = current_set.split();
        if(!definitely(split_sets.first.is_empty())) { _append_initial_set(working_sets,initial_time,split_sets.first); }
        if(!definitely(split_sets.second.is_empty())) { _append_initial_set(working_sets,initial_time,split_sets.second); }
    } else {
        working_sets.push_back(make_pair(initial_time,current_set));
    }
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

    // Track the previous step size used to properly find the starting step size
    StepSizeType previous_step_size = 0;

    ProgressIndicator initials_indicator(working_sets.size());
    ProgressIndicator time_indicator(maximum_time.get_d());

    while(!working_sets.empty()) {
        TimedEnclosureType current_timed_set=working_sets.back();
        working_sets.pop_back();
        TimeStepType current_time=current_timed_set.first;
        EnclosureType current_set_model=current_timed_set.second;
        FloatDPUpperBound current_set_radius=current_set_model.euclidean_set().bounding_box().radius();

        ARIADNE_LOG_PRINTLN("#w="<<std::setw(4)<<std::left<<working_sets.size()+1
                                 <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                                 <<" t="<<std::setw(7)<<std::fixed<<current_time.get_d()
                                 <<" p="<<std::setw(3)<<std::left<<current_set_model.number_of_parameters()
                                 <<" r="<<std::setw(6)<<current_set_model.radius()
                                 <<" ri="<<current_set_model.bounding_box().euclidean_set().radii()
                                 <<" c="<<current_set_model.centre());

        if(definitely(current_time>=maximum_time)) {
            final_sets.adjoin(current_set_model);
        } else if(semantics == Semantics::UPPER && configuration().enable_subdivisions() && decide(current_set_radius>configuration().maximum_enclosure_radius())) {
            // Subdivide
            List< EnclosureType > subdivisions=subdivide(current_set_model);
            for(SizeType i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_pair(current_time,subdivided_set_model));
            }
        } else if(semantics == Semantics::LOWER && configuration().enable_premature_termination() && decide(current_set_radius>configuration().maximum_enclosure_radius())) {
            ARIADNE_WARN("Terminating lower evolution at time " << current_time << " and set " << current_set_model << " due to maximum radius being exceeded.")
        } else {
            this->_evolution_step(working_sets,final_sets,reach_sets,intermediate_sets, current_timed_set,previous_step_size,maximum_time, semantics,reach);
        }

        initials_indicator.update_current(final_sets.size());
        time_indicator.update_current(current_time.get_d());
        if (initials_indicator.final_value() > 1) { ARIADNE_LOG_SCOPE_PRINTHOLD("[" << time_indicator.symbol() << "] " << initials_indicator.percentage() << "% of sets, " << time_indicator.percentage() << "% of current set"); }
        else ARIADNE_LOG_SCOPE_PRINTHOLD("[" << time_indicator.symbol() << "] " << time_indicator.percentage() << "%");
    }

    ARIADNE_LOG_PRINTLN("Finished evolution");
}


Void
VectorFieldEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& working_timed_set,
                StepSizeType& previous_step_size,
                const TimeType& maximum_time,
                Semantics semantics,
                Bool reach) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    EnclosureType current_set;
    TimeStepType current_time;
    ARIADNE_LOG_PRINTLN_AT(1,"working_timed_set_model = "<<working_timed_set);
    make_lpair(current_time,current_set)=working_timed_set;

    ARIADNE_LOG_PRINTLN("current_time = "<<current_time);
    ARIADNE_LOG_PRINTLN("current_set_model = "<<current_set);

    ARIADNE_LOG_PRINTLN("box = "<<current_set.bounding_box());
    ARIADNE_LOG_PRINTLN("radius = "<<current_set.euclidean_set().bounding_box().radius());

    /////////////// Main Evolution ////////////////////////////////

    // Push inputs
    runner()->push(TaskInput<VectorFieldEvolver>(system().dynamic_function(), current_set, current_time, previous_step_size));
    // Pull outputs
    auto out = runner()->pull();
    // Save outputs
    reach_sets.adjoin(out.reach);
    intermediate_sets.adjoin(out.evolve);
    working_sets.push_back(make_pair(out.time,out.evolve));
    previous_step_size = out.step_size_used;
}

Task<VectorFieldEvolver>::Task() : ParameterSearchTaskBase<C>("stp") { }

auto Task<VectorFieldEvolver>::run(TaskInput<C> const& in, Configuration<C> const& cfg) const -> TaskOutput<C> {
    LabelledEnclosure next_set = in.current_set;
    LabelledEnclosure reach_set = in.current_set;
    Dyadic next_time = in.current_time;
    Dyadic chosen_step_size = cfg.maximum_step_size();

    if(cfg.enable_reconditioning() && possibly(norm(next_set.state_function().errors()) > cfg.maximum_spacial_error())) {
        ARIADNE_LOG_PRINTLN_AT(1,"reconditioning from errors "<<next_set.state_function().errors());
        next_set.recondition();
    }

    auto set_bounds=cast_exact_box(next_set.euclidean_set().bounding_box());
    ARIADNE_LOG_PRINTLN_VAR_AT(1, set_bounds);

    FlowStepModelType flow_model = cfg.integrator().flow_step(in.dynamic, set_bounds, in.previous_step_size,chosen_step_size);

    ARIADNE_LOG_PRINTLN_VAR_AT(1, chosen_step_size);
    ARIADNE_LOG_PRINTLN_VAR_AT(1, flow_model);
    next_time += chosen_step_size;
    ARIADNE_LOG_PRINTLN_VAR_AT(1, next_time);
    reach_set.apply_full_reach_step(flow_model);
    ARIADNE_LOG_PRINTLN_VAR_AT(1, reach_set);
    next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
    ARIADNE_LOG_PRINTLN_VAR_AT(1, next_set);
    return TaskOutput<C>(next_set, reach_set, next_time, chosen_step_size);
}

}  // namespace Ariadne

