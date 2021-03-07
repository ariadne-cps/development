/***************************************************************************
 *            dynamics/inclusion_vector_field.hpp
 *
 *  Copyright  2018-20  Luca Geretti
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

/*! \file dynamics/inclusion_vector_field.hpp
 *  \brief Differential inclusion continuous dynamics system class.
 */

#ifndef ARIADNE_INCLUSION_VECTOR_FIELD_HPP
#define ARIADNE_INCLUSION_VECTOR_FIELD_HPP

#include "function/function.hpp"
#include "geometry/set_interface.hpp"
#include "geometry/grid.hpp"
#include "utility/identifier.hpp"
#include "symbolic/expression.decl.hpp"

namespace Ariadne {

class MissingInputException : public std::runtime_error {
  public:
    MissingInputException(const String& str) : std::runtime_error(str) { }
};

class UnusedInputException : public std::runtime_error {
  public:
    UnusedInputException(const String& str) : std::runtime_error(str) { }
};

class FunctionArgumentsMismatchException : public std::runtime_error {
  public:
    FunctionArgumentsMismatchException(const String& str) : std::runtime_error(str) { }
};

class Enclosure;
class InclusionVectorFieldEvolver;

//! \brief A vector field in Euclidean space, with noisy inputs identifying a family of trajectories.
class InclusionVectorField
{
  public:
    //! \brief The type used to represent time.
    typedef Real TimeType;
    //! \brief The type used to represent real numbers.
    typedef Real RealType;
    //! \brief The type used to describe the state space.
    typedef EuclideanSpace StateSpaceType;
    //! \brief The generic type used to compute the system evolution.
    typedef InclusionVectorFieldEvolver EvolverType;
    typedef Enclosure EnclosureType;
  public:
    //! \brief Construct from expressions with dotted \a dynamics and a given range for the \a inputs.
    InclusionVectorField(DottedRealAssignments const& dynamics, RealVariablesBox const& inputs);
    //! \brief Construct from a \a function for the dynamics and a given range for the \a inputs.
    InclusionVectorField(EffectiveVectorMultivariateFunction const& function, BoxDomainType const& inputs);
    virtual ~InclusionVectorField() = default;
    virtual InclusionVectorField* clone() const { return new InclusionVectorField(*this); }

    SizeType dimension() const { return _function.result_size(); }
    SizeType number_of_inputs() const { return _inputs.size(); }
    RealSpace state_space() const;
    Grid grid() const { return Grid(_function.result_size()); }

    //! \brief If the dynamics is affine in the inputs.
    Bool is_input_affine() const;
    //! \brief If the inputs are additive (includes the case of constant multipliers).
    Bool is_input_additive() const;

    const EffectiveVectorMultivariateFunction& function() const { return _function; }
    const BoxDomainType& inputs() const { return _inputs; }

    //! \brief Return the dynamics component obtained by setting the input radius to zero.
    EffectiveVectorMultivariateFunction noise_independent_component() const;
    //! \brief Return the dynamics components given by the derivatives for each input.
    Vector<EffectiveVectorMultivariateFunction> input_derivatives() const;

    friend OutputStream& operator<<(OutputStream& os, const InclusionVectorField& vf) {
        return os << "InclusionVectorField( " << vf.function() << ", " << vf.inputs() << " )"; }
  private:
    Void _transform_and_assign(EffectiveVectorMultivariateFunction const& function, BoxDomainType const& inputs);
  private:
    EffectiveVectorMultivariateFunction _function;
    BoxDomainType _inputs;
    List<Identifier> _variable_names;
};

} // namespace Ariadne

#endif // ARIADNE_INCLUSION_VECTOR_FIELD_HPP
