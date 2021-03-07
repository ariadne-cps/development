/***************************************************************************
 *            dynamics/vector_field_evolver.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/*! \file dynamics/vector_field_evolver.hpp
 *  \brief Evolver for vector_field systems.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_HPP
#define ARIADNE_VECTOR_FIELD_EVOLVER_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "utility/tuple.hpp"

#include "dynamics/vector_field.hpp"
#include "function/function_interface.hpp"
#include "configuration/configuration_interface.hpp"
#include "solvers/integrator_interface.hpp"
#include "solvers/integrator.hpp"
#include "dynamics/evolver_base.hpp"
#include "output/logging.hpp"
#include "concurrency/task_ranking_space.hpp"
#include "configuration/searchable_configuration.hpp"
#include "concurrency/task_runner_interface.hpp"
#include "concurrency/task_runner.tpl.hpp"
#include "configuration/configuration_property.hpp"

namespace Ariadne {

class VectorField;
template<class ES> class Orbit;

class VectorFieldEvolver;
template<> struct Configuration<VectorFieldEvolver>;

//! \brief A class for computing the evolution of a vector_field system.
//!
//! The actual evolution steps are performed by the Integrator class.
class VectorFieldEvolver
    : public EvolverBase<VectorField,LabelledEnclosure,typename VectorField::TimeType>,
      public TaskRunnable<VectorFieldEvolver>
{
  public:
    typedef VectorField SystemType;
    typedef Configuration<VectorFieldEvolver> ConfigurationType;
    typedef typename VectorField::TimeType TimeType;
    typedef Dyadic TimeStepType;
    typedef TimeType TerminationType;
    typedef LabelledEnclosure EnclosureType;
    typedef Pair<TimeStepType, EnclosureType> TimedEnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
  public:

    //! \brief Construct from parameters and an integrator to compute the flow.
    VectorFieldEvolver(SystemType const& system, ConfigurationType const& configuration);

    //! \brief Make a dynamically-allocated copy.
    VectorFieldEvolver* clone() const override;

    //! \brief Get the internal system.
    SystemType const& system() const override { return *_sys_ptr; }

    //! \brief Make an enclosure from a user set.
    EnclosureType enclosure(RealBox const&) const;
    EnclosureType enclosure(RealBox const&, EnclosureConfiguration const&) const;

    //! \brief Make an enclosure from a user set with variables.
    EnclosureType enclosure(RealVariablesBox const&) const;
    EnclosureType enclosure(RealVariablesBox const&, EnclosureConfiguration const&) const;

    //! \brief Make an enclosure from a computed box set.
    EnclosureType enclosure(ExactBoxType const&) const;
    EnclosureType enclosure(ExactBoxType const&, EnclosureConfiguration const&) const;

    //!@}

    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_set, const TimeType& time, Semantics semantics=Semantics::UPPER) const override;

    using EvolverBase< VectorField, EnclosureType, TerminationType >::evolve;
    using EvolverBase< VectorField, EnclosureType, TerminationType >::reach;

    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,Semantics::UPPER,false);
        return final; }

    //! \brief Compute an approximation to the reachable set under upper semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,Semantics::UPPER,true);
        return reachable; }
    //!@}

protected:

    virtual Void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, Bool reach) const override;

    virtual Void _evolution_step(List< TimedEnclosureType >& working_sets,
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                 const TimedEnclosureType& current_set, StepSizeType& last_step_size, const TimeType& time,
                                 Semantics semantics, Bool reach) const;

    virtual Void _append_initial_set(List<TimedEnclosureType>& working_sets, const TimeStepType& initial_time, const EnclosureType& current_set) const;

  private:
    SharedPointer<SystemType> _sys_ptr;
};

//! \brief Configuration for a VectorFieldEvolver, essentially for controlling the accuracy of continuous evolution methods.
template<> struct Configuration<VectorFieldEvolver> final : public SearchableConfiguration {
    typedef Configuration<VectorFieldEvolver> C;
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<IntegratorInterface> IntegratorProperty;

    Configuration() {
        add_property("enable_premature_termination",BooleanConfigurationProperty(false));
        add_property("enable_reconditioning",BooleanConfigurationProperty(false));
        add_property("enable_subdivisions",BooleanConfigurationProperty(false));
        add_property("integrator", InterfaceListConfigurationProperty<IntegratorInterface>(TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator>())));
        add_property("maximum_enclosure_radius",RealTypeProperty(ExactDouble::infinity(),Log10SearchSpaceConverter<RealType>()));
        add_property("maximum_spacial_error",RealTypeProperty(ExactDouble::infinity(),Log10SearchSpaceConverter<RealType>()));
        add_property("maximum_step_size",RealTypeProperty(ExactDouble::infinity(),Log2SearchSpaceConverter<RealType>()));
    }

    //! \brief Enable premature termination of lower evolution
    Bool const& enable_premature_termination() const { return at<BooleanConfigurationProperty>("enable_premature_termination").get(); }
    C& set_enable_premature_termination(Bool const& value) { at<BooleanConfigurationProperty>("enable_premature_termination").set(value); return *this; }

    //! \brief Enable reconditioning of basic sets
    Bool const& enable_reconditioning() const { return at<BooleanConfigurationProperty>("enable_reconditioning").get(); }
    C& set_enable_reconditioning(Bool const& value) { at<BooleanConfigurationProperty>("enable_reconditioning").set(value); return *this; }
    C& set_both_enable_reconditioning() { at<BooleanConfigurationProperty>("enable_reconditioning").set_both(); return *this; }

    //! \brief Enable subdivisions for upper evolution
    Bool const& enable_subdivisions() const { return at<BooleanConfigurationProperty>("enable_subdivisions").get(); }
    C& set_enable_subdivisions(Bool const& value) { at<BooleanConfigurationProperty>("enable_subdivisions").set(value); return *this; }

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType const& maximum_step_size() const { return at<RealTypeProperty>("maximum_step_size").get(); }
    C& set_maximum_step_size(ApproximateRealType const& value) { at<RealTypeProperty>("maximum_step_size").set(cast_exact(value)); return *this; }
    C& set_maximum_step_size(ApproximateRealType const& lower, ApproximateRealType const& upper) { at<RealTypeProperty>("maximum_step_size").set(cast_exact(lower),cast_exact(upper)); return *this; }

    //! \brief The maximum allowable approximation error for reconditioning
    RealType const& maximum_spacial_error() const { return at<RealTypeProperty>("maximum_spacial_error").get(); }
    C& set_maximum_spacial_error(ApproximateRealType const& value) { at<RealTypeProperty>("maximum_spacial_error").set(cast_exact(value)); return *this; }
    C& set_maximum_spacial_error(ApproximateRealType const& lower, ApproximateRealType const& upper) { at<RealTypeProperty>("maximum_spacial_error").set(cast_exact(lower),cast_exact(upper)); return *this; }

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType const& maximum_enclosure_radius() const { return at<RealTypeProperty>("maximum_enclosure_radius").get(); }
    C& set_maximum_enclosure_radius(ApproximateRealType const& value) { at<RealTypeProperty>("maximum_enclosure_radius").set(cast_exact(value)); return *this; }

    //! \brief The integrator to be used.
    IntegratorInterface const& integrator() const { return at<IntegratorProperty>("integrator").get(); }
    C& set_integrator(IntegratorInterface const& integrator) { at<IntegratorProperty>("integrator").set(integrator); return *this; }
    C& set_integrator(SharedPointer<IntegratorInterface> const& integrator) { at<IntegratorProperty>("integrator").set(integrator); return *this; }
};

template<> struct TaskInput<VectorFieldEvolver> {
    TaskInput(EffectiveVectorMultivariateFunction const& dynamic_, LabelledEnclosure const& current_set_,
              Dyadic const& current_time_, Dyadic const& previous_step_size_) :
            dynamic(dynamic_), current_set(current_set_), current_time(current_time_), previous_step_size(previous_step_size_) { }
    EffectiveVectorMultivariateFunction const& dynamic;
    LabelledEnclosure const& current_set;
    Dyadic const& current_time;
    Dyadic const& previous_step_size;
};

template<> struct TaskOutput<VectorFieldEvolver> {
    TaskOutput(LabelledEnclosure const& evolve_, LabelledEnclosure const& reach_, Dyadic const& time_, Dyadic const& step_size_used_) :
            evolve(evolve_), reach(reach_), time(time_), step_size_used(step_size_used_) { }
    LabelledEnclosure const evolve;
    LabelledEnclosure const reach;
    Dyadic const time;
    Dyadic const step_size_used;
};

template<> struct TaskObjective<VectorFieldEvolver> {
    TaskObjective(RealVariable const& v, PositiveFloatDPUpperBound const& r, Dyadic const& t) : variable(v), radius(r), time(t) { }
    RealVariable variable;
    PositiveFloatDPUpperBound const radius;
    Dyadic const time;
};

template<> struct Task<VectorFieldEvolver> final: public ParameterSearchTaskBase<VectorFieldEvolver> {
    typedef VectorFieldEvolver C;
    Task();
    TaskOutput<C> run(TaskInput<C> const& in, Configuration<C> const& cfg) const override;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_HPP
