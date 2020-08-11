/***************************************************************************
 *            function/function_interface.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file function/function_interface.hpp
 *  \brief Interface for functions for which derivatives can be computed.
 */

#ifndef ARIADNE_FUNCTION_INTERFACE_HPP
#define ARIADNE_FUNCTION_INTERFACE_HPP

#include <iosfwd>

#include "utility/declarations.hpp"
#include "function/function.decl.hpp"

namespace Ariadne {

static const Int SMOOTH=255;

template<class S> struct ElementTraits;
template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;
template<class S> using ElementKind = typename ElementTraits<S>::Kind;
template<class S> using ElementSizeType = decltype(declval<S>().dimension());
template<class S> using ElementIndexType = decltype(declval<S>().dimension());

template<class... ARGS> struct DomainOfTypedef;
template<> struct DomainOfTypedef<RealScalar> { typedef IntervalDomainType Type; };
template<> struct DomainOfTypedef<RealVector> { typedef BoxDomainType Type; };
template<class... ARGS> using DomainOfType = typename DomainOfTypedef<ARGS...>::Type;


template<class SIG> struct SignatureTraits;
template<> struct SignatureTraits<Real(Real)> {
    typedef Real ResultKind; typedef Real ArgumentKind;
    typedef RealDomain DomainType; typedef RealDomain CodomainType;
    typedef IntervalDomainType BoundedDomainType; typedef IntervalDomainType BoundedCodomainType;
    typedef IntervalRangeType BoundedRangeType;
};
template<> struct SignatureTraits<Real(RealVector)> {
    typedef Real ResultKind; typedef RealVector ArgumentKind;
    typedef EuclideanDomain DomainType; typedef RealDomain CodomainType;
    typedef BoxDomainType BoundedDomainType; typedef IntervalDomainType BoundedCodomainType;
    typedef IntervalRangeType BoundedRangeType;
};
template<> struct SignatureTraits<RealVector(Real)> {
    typedef RealVector ResultKind; typedef Real ArgumentKind;
    typedef RealDomain DomainType; typedef EuclideanDomain CodomainType;
    typedef IntervalDomainType BoundedDomainType; typedef BoxDomainType BoundedCodomainType;
    typedef BoxRangeType BoundedRangeType;
};
template<> struct SignatureTraits<RealVector(RealVector)> {
    typedef RealVector ResultKind; typedef RealVector ArgumentKind;
    typedef EuclideanDomain DomainType; typedef EuclideanDomain CodomainType;
    typedef BoxDomainType BoundedDomainType; typedef BoxDomainType BoundedCodomainType;
    typedef BoxRangeType BoundedRangeType;
};


template<class P, class SIG> class FunctionInterface;

template<class P, class... ARGS> using ScalarFunctionInterface = FunctionInterface<P,RealScalar(ARGS...)>;
template<class P, class... ARGS> using VectorFunctionInterface = FunctionInterface<P,RealVector(ARGS...)>;
template<class P> using ScalarUnivariateFunctionInterface = ScalarFunctionInterface<P,RealScalar>;
template<class P> using VectorUnivariateFunctionInterface = VectorFunctionInterface<P,RealScalar>;
template<class P> using ScalarMultivariateFunctionInterface = ScalarFunctionInterface<P,RealVector>;
template<class P> using VectorMultivariateFunctionInterface = VectorFunctionInterface<P,RealVector>;



template<class P,class... ARGS>
class VectorOfFunctionInterface
{
};

template<class... ARGS> class VectorOfFunctionInterface<ApproximateTag,ARGS...> {
  public:
    virtual ~VectorOfFunctionInterface<ApproximateTag,ARGS...>() = default;
    virtual ScalarFunctionInterface<ApproximateTag,ARGS...>* _get(SizeType i) const = 0;
};

template<class... ARGS> class VectorOfFunctionInterface<ValidatedTag,ARGS...>
    : public virtual VectorOfFunctionInterface<ApproximateTag,ARGS...>
{
  public:
    virtual ScalarFunctionInterface<ValidatedTag,ARGS...>* _get(SizeType i) const override = 0;
};

template<class... ARGS> class VectorOfFunctionInterface<EffectiveTag,ARGS...>
    : public virtual VectorOfFunctionInterface<ValidatedTag,ARGS...>
{
  public:
    virtual ScalarFunctionInterface<EffectiveTag,ARGS...>* _get(SizeType i) const override = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for functions \f$\R^n\rightarrow\R^m\f$ which can be evaluated.
template<class P, class SIG> class FunctionInterface;


template<class SIG>
class FunctionInterface<Void,SIG>
{
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<DomainType> ArgumentSizeType;
    typedef ElementSizeType<CodomainType> ResultSizeType;
    typedef ElementIndexType<DomainType> ArgumentIndexType;

    virtual ~FunctionInterface() = default;
    virtual ArgumentSizeType argument_size() const = 0;
    virtual ResultSizeType result_size() const = 0;
//    virtual DomainType const domain() const = 0;
//    virtual CodomainType const codomain() const = 0;

    virtual OutputStream& repr(OutputStream& os) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
  public:
    virtual FunctionInterface<Void,SIG>* _clone() const = 0;
    inline FunctionInterface<Void,SIG>* _copy() const { return this->_clone(); }
  public:
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionInterface<Void,SIG>& f) { return f._write(os); }
};

//! \ingroup FunctionModule
//! \brief Interface for functions \f$\R^n\rightarrow\R^m\f$ which can only be evaluated approximately.
//! \sa \ref FunctionInterface "FunctionInterface<P,SIG>"
template<class SIG>
class FunctionInterface<ApproximateTag,SIG>
    : public virtual FunctionInterface<Void,SIG>
{
    using P=ApproximateTag;
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
  public:
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    virtual Result<FloatDPApproximation> _evaluate(const Argument<FloatDPApproximation>& x) const = 0;
    virtual Result<FloatMPApproximation> _evaluate(const Argument<FloatMPApproximation>& x) const = 0;
    virtual Result<Differential<FloatDPApproximation>> _evaluate(const Argument< Differential<FloatDPApproximation> >& x) const = 0;
    virtual Result<Differential<FloatMPApproximation>> _evaluate(const Argument< Differential<FloatMPApproximation> >& x) const = 0;
    virtual Result<TaylorModel<ApproximateTag,FloatDP>> _evaluate(const Argument< TaylorModel<ApproximateTag,FloatDP> >& x) const = 0;
    virtual Result<TaylorModel<ApproximateTag,FloatMP>> _evaluate(const Argument< TaylorModel<ApproximateTag,FloatMP> >& x) const = 0;
    virtual Result<Formula<ApproximateNumber>> _evaluate(const Argument< Formula<ApproximateNumber> >& x) const = 0;
    virtual Result<ElementaryAlgebra<ApproximateNumber>> _evaluate(const Argument< ElementaryAlgebra<ApproximateNumber> >& x) const = 0;

    virtual FunctionInterface<P,SIG>* _clone() const = 0;
    inline FunctionInterface<P,SIG>* _copy() const { return this->_clone(); }
    virtual FunctionInterface<P,SIG>* _derivative(ElementIndexType<D> i) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for functions \f$\R^n\rightarrow\R\f$ which can be evaluated with guaranteed bounds on the error.
//! \sa \ref FunctionInterface "FunctionInterface<P,SIG>"
template<class SIG>
class FunctionInterface<ValidatedTag,SIG>
    : public virtual FunctionInterface<ApproximateTag,SIG>
{
    using P = ValidatedTag;
    using AP = ApproximateTag;
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
  public:
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionInterface<AP,SIG>::_evaluate;
    virtual Result<FloatDPBounds> _evaluate(const Argument<FloatDPBounds>& x) const = 0;
    virtual Result<FloatMPBounds> _evaluate(const Argument<FloatMPBounds>& x) const = 0;
    virtual Result<Differential<FloatDPBounds>> _evaluate(const Argument< Differential<FloatDPBounds> >& x) const = 0;
    virtual Result<Differential<FloatMPBounds>> _evaluate(const Argument< Differential<FloatMPBounds> >& x) const = 0;
    virtual Result<TaylorModel<ValidatedTag,FloatDP>> _evaluate(const Argument< TaylorModel<ValidatedTag,FloatDP> >& x) const = 0;
    virtual Result<TaylorModel<ValidatedTag,FloatMP>> _evaluate(const Argument< TaylorModel<ValidatedTag,FloatMP> >& x) const = 0;
    virtual Result<TaylorModel<ValidatedTag,FloatDPUpperInterval>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatDPUpperInterval>>& x) const = 0;
    virtual Result<TaylorModel<ValidatedTag,FloatMPUpperInterval>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatMPUpperInterval>>& x) const = 0;

    virtual Result<Formula<ValidatedNumber>> _evaluate(const Argument< Formula<ValidatedNumber> >& x) const = 0;
    virtual Result<ElementaryAlgebra<ValidatedNumber>> _evaluate(const Argument< ElementaryAlgebra<ValidatedNumber> >& x) const = 0;

    virtual Result<ScalarMultivariateFunction<ValidatedTag>> _evaluate(const Argument< ScalarMultivariateFunction<ValidatedTag> >& x) const = 0;

    inline Result<FloatDPBounds> _evaluate(const Argument<FloatDPValue>& x) const {
        return this->_evaluate(Argument<FloatDPBounds>(x)); }
    inline Result<FloatMPBounds> _evaluate(const Argument<FloatMPValue>& x) const {
        return this->_evaluate(Argument<FloatMPBounds>(x)); }

    virtual FunctionInterface<P,SIG>* _clone() const = 0;
    inline FunctionInterface<P,SIG>* _copy() const { return this->_clone(); }
    virtual FunctionInterface<P,SIG>* _derivative(ElementIndexType<D> i) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref FunctionInterface "FunctionInterface<P,SIG>"
template<class SIG>
class FunctionInterface<EffectiveTag,SIG>
    : public virtual FunctionInterface<ValidatedTag,SIG>
{
    using WP = ValidatedTag;
    using P = EffectiveTag;
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
  public:
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionInterface<WP,SIG>::_evaluate;
    virtual Result<Real> _evaluate(const Argument<Real>& x) const = 0;
    virtual Result<ElementaryAlgebra<Real>> _evaluate(const Argument<ElementaryAlgebra<Real>>& x) const = 0;
    virtual Result<Formula<Real>> _evaluate(const Argument<Formula<Real>>& x) const = 0;
    virtual Result<ElementaryAlgebra<EffectiveNumber>> _evaluate(const Argument<ElementaryAlgebra<EffectiveNumber>>& x) const = 0;
    virtual Result<Formula<EffectiveNumber>> _evaluate(const Argument<Formula<EffectiveNumber>>& x) const = 0;

    virtual FunctionInterface<P,SIG>* _clone() const = 0;
    inline FunctionInterface<P,SIG>* _copy() const { return this->_clone(); }
    virtual FunctionInterface<P,SIG>* _derivative(ElementIndexType<D> i) const = 0;
};


template<class I> class VectorInterface {
    ~VectorInterface() = default;
    virtual I* _get(SizeType i) const = 0;
};

template<class P> class FunctionFactoryInterface;

template<> class FunctionFactoryInterface<ValidatedTag>
{
    using P = ValidatedTag;
    typedef BoxDomainType DomainType;
  public:
    virtual FunctionFactoryInterface<P>* clone() const = 0;
    inline FunctionFactoryInterface<P>* _copy() const { return this->clone(); }
    virtual OutputStream& _write(OutputStream& os) const = 0;
    inline ScalarMultivariateFunction<P> create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const;
    inline VectorMultivariateFunction<P> create(const BoxDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const;
    inline ScalarMultivariateFunction<P> create_zero(const BoxDomainType& domain) const;
    inline VectorMultivariateFunction<P> create_identity(const BoxDomainType& domain) const;
  private:
    virtual ScalarMultivariateFunctionInterface<P>* _create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const = 0;
    virtual VectorMultivariateFunctionInterface<P>* _create(const BoxDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const = 0;
  public:
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionFactoryInterface<P>& factory) {
        return factory._write(os); }
  public:
    virtual ~FunctionFactoryInterface<P>() = default;
};

} // namespace Ariadne

#endif
