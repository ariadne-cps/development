/***************************************************************************
 *            inclusion_evolver.hpp
 *
 *  Copyright  2008-20  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

/*! \file inclusion_evolver.hpp
 *  \brief Evolver for differential inclusion dynamics.
 */

#ifndef ARIADNE_INCLUSION_EVOLVER_HPP
#define ARIADNE_INCLUSION_EVOLVER_HPP

#include "utility/typedefs.hpp"
#include "utility/attribute.hpp"
#include "numeric/numeric.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/algebra.hpp"
#include "function/domain.hpp"
#include "function/function_patch.hpp"
#include "function/function_model.hpp"
#include "function/formula.hpp"
#include "function/symbolic_function.hpp"
#include "symbolic/expression_set.hpp"
#include "output/logging.hpp"
#include "solvers/integrator_interface.hpp"
#include "solvers/inclusion_integrator.hpp"
#include "solvers/configuration_interface.hpp"
#include "inclusion_vector_field.hpp"

namespace Ariadne {

class Real;

struct StepSize : public Attribute<StepSizeType> { };
struct NumberOfStepsBetweenSimplifications : public Attribute<Nat> { };
struct NumberOfVariablesToKeep : public Attribute<Nat> { };

Generator<StepSize> step_size;
Generator<NumberOfStepsBetweenSimplifications> number_of_steps_between_simplifications;
Generator<NumberOfVariablesToKeep> number_of_variables_to_keep;

using ThresholdSweeperDP = ThresholdSweeper<FloatDP>;
using GradedSweeperDP = GradedSweeper<FloatDP>;
using GradedThresholdSweeperDP = GradedThresholdSweeper<FloatDP>;
using SweeperDP = Sweeper<FloatDP>;

using TimeStepType = Dyadic;

BoxDomainType initial_ranges_to_box(RealVariablesBox const& var_ranges);

FloatDP volume(Vector<ApproximateIntervalType> const& box);

class InclusionEvolverState;

class ReconditionerInterface {
  public:
    virtual Bool must_reduce_parameters(InclusionEvolverState const& state) const = 0;
    virtual Bool must_incorporate_errors(InclusionEvolverState const& state) const = 0;
    virtual Void reduce_parameters(ValidatedVectorMultivariateFunctionPatch& phi) const = 0;
    virtual ValidatedVectorMultivariateFunctionPatch incorporate_errors(ValidatedVectorMultivariateFunctionPatch const& Phi) const = 0;
    virtual Void update_from(InclusionEvolverState const& state) = 0;
  public:
    virtual ReconditionerInterface* clone() const = 0;
    virtual ~ReconditionerInterface() = default;
};

class LohnerReconditioner : public ReconditionerInterface {
    Nat _number_of_variables;
    Nat _number_of_inputs;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_parameters_to_keep;
    ExactDouble _ratio_of_parameters_to_keep;
public:
    LohnerReconditioner(Nat number_of_variables, Nat number_of_inputs, Nat number_of_steps_between_simplifications_, ApproximateDouble ratio_of_parameters_to_keep)
        : _number_of_variables(number_of_variables),
          _number_of_inputs(number_of_inputs),
          _number_of_steps_between_simplifications(number_of_steps_between_simplifications_),
          _number_of_parameters_to_keep(USHRT_MAX),
          _ratio_of_parameters_to_keep(cast_exact(ratio_of_parameters_to_keep)) { }
    virtual LohnerReconditioner* clone() const override { return new LohnerReconditioner(*this); }
    virtual ValidatedVectorMultivariateFunctionPatch incorporate_errors(ValidatedVectorMultivariateFunctionPatch const& f) const override;
    virtual Void reduce_parameters(ValidatedVectorMultivariateFunctionPatch& f) const override;
    virtual Bool must_reduce_parameters(InclusionEvolverState const& state) const override;
    virtual Bool must_incorporate_errors(InclusionEvolverState const& state) const override;
    virtual Void update_from(InclusionEvolverState const& state) override;
};

class ReconditionerHandle {
private:
    SharedPointer<ReconditionerInterface> _impl;
public:
    ReconditionerHandle(ReconditionerInterface const& other) : _impl(other.clone()) { }
    ReconditionerHandle(ReconditionerHandle const& other) : _impl(other._impl) { }

    ReconditionerHandle& operator=(ReconditionerHandle const& other) { _impl = other._impl; return *this; }

    Bool must_reduce_parameters(InclusionEvolverState const& state) const { return _impl->must_reduce_parameters(state); }
    Bool must_incorporate_errors(InclusionEvolverState const& state) const { return _impl->must_incorporate_errors(state); }

    Void reduce_parameters(ValidatedVectorMultivariateFunctionPatch& phi) const { _impl->reduce_parameters(phi); }
    ValidatedVectorMultivariateFunctionPatch incorporate_errors(ValidatedVectorMultivariateFunctionPatch const& Phi) const { return _impl->incorporate_errors(Phi); }

    Void update_from(InclusionEvolverState const& state) {_impl->update_from(state); }

    virtual ~ReconditionerHandle() = default;
};

class InclusionEvolverConfiguration;

class InclusionEvolver {
  public:
    typedef InclusionEvolverConfiguration ConfigurationType;
    typedef InclusionVectorField SystemType;
  protected:
    SystemType _system;
    SweeperDP _sweeper;
    SharedPointer<IntegratorInterface> _integrator;
    ReconditionerHandle _reconditioner;
    SharedPointer<ConfigurationType> _configuration;
  public:
    InclusionEvolver(SystemType const& system, SweeperDP const& sweeper, IntegratorInterface const& integrator, ReconditionerHandle const& reconditioner);

    //!@{
    //! \name Configuration for the class.
    //! \brief A reference to the configuration controlling the evolution.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

  public:
    List<ValidatedVectorMultivariateFunctionPatch> reach(BoxDomainType const& initial, Real const& T);
  private:
    Void _recondition_and_update(ValidatedVectorMultivariateFunctionPatch& function, InclusionEvolverState& state);
};

//! \brief Configuration for an InclusionEvolver, essentially for controlling the accuracy of continuous evolution methods.
class InclusionEvolverConfiguration : public ConfigurationInterface
{
  public:
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;

    //! \brief Default constructor gives reasonable values.
    InclusionEvolverConfiguration();

    virtual ~InclusionEvolverConfiguration() = default;

  private:

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType _maximum_step_size;

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_enclosure_radius;

    //! \brief Enable reduction of the parameters (true by default).
    Bool _enable_parameter_reduction;

    //! \brief The approximations allowed in calculating the flow for each step.
    //! At least one approximation is required.
    List<InputApproximation> _approximations;

  public:

    const RealType& maximum_step_size() const { return _maximum_step_size; }
    Void maximum_step_size(const ApproximateRealType value) { _maximum_step_size = cast_exact(value); }

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    Void maximum_enclosure_radius(const ApproximateRealType value) { _maximum_enclosure_radius = cast_exact(value); }

    const Bool& enable_parameter_reduction() const { return _enable_parameter_reduction; }
    Void enable_parameter_reduction(const Bool value) { _enable_parameter_reduction = value; }

    List<InputApproximation> const& approximations() const { return _approximations; }
    Void approximations(List<InputApproximation> const& value) { assert(value.size()>0); _approximations = value; }

  public:

    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne;

#endif // ARIADNE_INCLUSION_EVOLVER_HPP
