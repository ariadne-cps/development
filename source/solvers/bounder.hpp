/***************************************************************************
 *            solvers/bounder.hpp
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

/*! \file solvers/bounder.hpp
 *  \brief Classes for finding the bounds of the flow for differential equations.
 */

#ifndef ARIADNE_BOUNDER_HPP
#define ARIADNE_BOUNDER_HPP

#include "utility/typedefs.hpp"
#include "utility/attribute.hpp"
#include "algebra/algebra.hpp"
#include "function/domain.hpp"
#include "function/function_model.hpp"
#include "output/logging.hpp"
#include "integrator_interface.hpp"

namespace Ariadne {

class BoundingNotFoundException : public std::runtime_error {
  public:
    BoundingNotFoundException(const String& str) : std::runtime_error(str) { }
};

//! \ingroup SolverModule EvaluationModule
//! \brief Interface for classes calculating the bounds of a flow.
class BounderInterface {
  public:
    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, StepSizeType const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& vector_field, BoxDomainType const& state_domain, BoxDomainType const& parameter_domain, StepSizeType const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, BoxDomainType const& parameter_domain, StepSizeType const& suggested_time_step) const = 0;

  public:
    virtual Void _write(OutputStream& os) const = 0;
    virtual BounderInterface* clone() const = 0;

    friend inline OutputStream& operator<<(OutputStream& os, BounderInterface const& bounder) {
        bounder._write(os); return os; }

    virtual ~BounderInterface() = default;
};

class BounderBase : public BounderInterface {
  public:
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const override = 0;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const override = 0;
};

class EulerBounder final : public BounderBase {
  public:
    using BounderBase::compute;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const override;
    virtual Void _write(OutputStream& os) const override { os << "EulerBounder"; }
    virtual EulerBounder* clone() const override { return new EulerBounder(*this); }
  private:
    // Compute the bounds on the reach set of dx/dt = f(x,t,a) starting at time t, for a suggested step size of h.
    Pair<StepSizeType,UpperBoxType> _compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const;
    // Compute a bounding box for the flow of f starting in D at T.lower() to time T.upper() with parameters A, assuming flow remains in B,
    // but widening D by INITIAL_WIDENING and the flow vectors by VECTOR_WIDENING.
    UpperBoxType _formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B, PositiveFloatDPValue INITIAL_WIDENING, PositiveFloatDPValue VECTOR_WIDENING) const;
    UpperBoxType _refinement(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const;
};

class BounderHandle {
  private:
    SharedPointer<BounderInterface> _impl;
  public:
    BounderHandle(BounderInterface const& bounder) : _impl(bounder.clone()) { }
    BounderHandle(BounderHandle const& other) : _impl(other._impl) { }
    BounderHandle& operator=(BounderHandle const& other) { _impl = other._impl; return *this; }

    operator BounderInterface const& () const { return *_impl; }

    Void _write(OutputStream& os) const { _impl->_write(os); }
    BounderHandle* clone() const { return new BounderHandle(*this); }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const {
        return _impl->compute(f,D,hsug);
    }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const {
        return _impl->compute(f,D,A,hsug);
    }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const {
        return _impl->compute(f,D,t,hsug);
    }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const {
        return _impl->compute(f,D,t,A,hsug);
    }

};

} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
