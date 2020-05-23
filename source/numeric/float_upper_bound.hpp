/***************************************************************************
 *            numeric/float_upper_bound.hpp
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

/*! \file numeric/float_upper_bound.hpp
 *  \brief Floating-point upper bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_UPPER_BOUND_HPP
#define ARIADNE_FLOAT_UPPER_BOUND_HPP

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
//! \sa UpperReal, FloatDP, FloatMP, Bounds, LowerBound.
template<class F> class UpperBound
    : public DefineDirectedGroupOperators<UpperBound<F>,LowerBound<F>>
    , public DefineDirectedGroupOperators<LowerBound<F>,UpperBound<F>>
    , public DefineDirectedComparisonOperators<UpperBound<F>,LowerBound<F>,LessTrait<UpperBound<F>>,EqualsTrait<UpperBound<F>>>
    , public DefineDirectedComparisonOperators<LowerBound<F>,UpperBound<F>,LessTrait<LowerBound<F>>,EqualsTrait<LowerBound<F>>>
    , public DefineConcreteGenericOperators<UpperBound<F>>
    , public DeclareFloatOperations<Approximation<F>>
{
  protected:
    typedef UpperTag P; typedef typename F::RoundingModeType RND; typedef typename F::PrecisionType PR;
  public:
    //! <p/>
    typedef UpperTag Paradigm;
    //! <p/>
    typedef UpperBound<F> NumericType;
    //! <p/>
    typedef ValidatedUpperNumber GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR PropertiesType;
  public:
    //! A upper bound of zero with precision \a pr.
    explicit UpperBound<F>(PrecisionType pr) : _u(0.0_x,pr) { }
    //! A upper bound with value \a u.
    explicit UpperBound<F>(RawType const& u) : _u(u) { }

    template<BuiltinIntegral N> UpperBound<F>(N n, PR pr) : UpperBound<F>(ExactDouble(n),pr) { }
    UpperBound<F>(const ExactDouble& d, PR pr) : _u(d,pr) { }
        UpperBound<F>(const TwoExp& t, PR pr) : _u(t,pr) { }
        UpperBound<F>(const Integer& z, PR pr) : _u(z,up,pr) { }
        UpperBound<F>(const Dyadic& w, PR pr) : _u(w,up,pr) { }
        UpperBound<F>(const Decimal& d, PR pr) : _u(d,up,pr) { }
        UpperBound<F>(const Rational& q, PR pr) : _u(q,up,pr) { }
        UpperBound<F>(const Real& r, PR pr);
        UpperBound<F>(const Value<F>& x, PR pr); // FIXME: Should not be necessary
        UpperBound<F>(const Bounds<F>& x, PR pr);
    UpperBound<F>(const UpperBound<F>& x, PR pr);
    //! A upper bound of type \p F from a generic upper bound \a y.
    UpperBound<F>(const ValidatedUpperNumber& y, PR pr);
    template<class FF> requires Constructible<F,FF,RND,PR>
        UpperBound<F>(const UpperBound<FF>& x, PR pr) : _u(x.raw(),up,pr) { }

    //! Convert from upper \em and lower bounds on a number.
    UpperBound<F>(Bounds<F> const& x);
    template<class FE> UpperBound<F>(Ball<F,FE> const& x);
    UpperBound<F>(Value<F> const& x);
    UpperBound<F>(Error<F> const& x); // FIXME: Remove

        UpperBound<F>& operator=(const Value<F>& x) { return *this=UpperBound<F>(x); }
    //! Assign from the upper bound \a y, keeping the same precision.
    UpperBound<F>& operator=(const ValidatedUpperNumber& y) { return *this=UpperBound<F>(y,this->precision()); }
    //! Create a upper bound from the generic upper bound \a y with the same precision as \a this.
    UpperBound<F> create(const ValidatedUpperNumber& y) const { return UpperBound<F>(y,this->precision()); }
    //! Create a lower bound from the generic lower bound \a y with the same precision as \a this.
    LowerBound<F> create(const ValidatedLowerNumber& y) const { return LowerBound<F>(y,this->precision()); }

    //! Downcast to a generic upper bound.
    operator ValidatedUpperNumber () const;

    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _u.precision(); }
    //! The compuational properties needed to create the upper bound; equivalent to the precision.
    PropertiesType properties() const { return _u.precision(); }
    //! Downcast to a generic upper bound.
    GenericType generic() const { return this->operator GenericType(); }
    //! The raw data used to represent the upper bound.
    RawType const& raw() const { return _u; }
    //! A mutable reference to the raw data used to represent the upper bound.
    RawType& raw() { return _u; }
    //! Over-approximate by a builtin double-precision value. DEPRECATED
    double get_d() const { return _u.get_d(); }
  public:
    friend Bool is_nan(UpperBound<F> const& x) {
        return is_nan(x._u); }

    friend UpperBound<F> max(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(max(x1._u,x2._u)); }
    friend UpperBound<F> min(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(min(x1._u,x2._u)); }
    friend Approximation<F> abs(UpperBound<F> const& x) {
        return abs(Approximation<F>(x)); }

    friend UpperBound<F> nul(UpperBound<F> const& x) {
        return UpperBound<F>(pos(x._u)); }
    friend UpperBound<F> pos(UpperBound<F> const& x) {
        return UpperBound<F>(pos(x._u)); }
    friend LowerBound<F> neg(UpperBound<F> const& x) {
        return LowerBound<F>(neg(x._u)); }
    friend UpperBound<F> hlf(UpperBound<F> const& x) {
        return UpperBound<F>(hlf(x._u)); }

    friend UpperBound<F> add(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(add(up,x1._u,x2._u)); }
    friend UpperBound<F> sub(UpperBound<F> const& x1, LowerBound<F> const& x2) {
        return UpperBound<F>(sub(up,x1._u,x2._l)); }

    friend PositiveUpperBound<F> sqrt(PositiveUpperBound<F> const& x);
    friend PositiveUpperBound<F> exp(UpperBound<F> const& x);
    friend UpperBound<F> log(PositiveUpperBound<F> const& x);
    friend UpperBound<F> atan(UpperBound<F> const& x);
    friend PositiveUpperBound<F> atan(PositiveUpperBound<F> const& x);

    friend ValidatedNegatedSierpinskian eq(UpperBound<F> const& x1, LowerBound<F> const& x2) {
        if(x1._u<x2._l) { return false; }
        else { return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE); } }
    friend ValidatedLowerKleenean lt(UpperBound<F> const& x1, LowerBound<F> const& x2) {
        if(x1._u< x2._l) { return true; }
        else { return ValidatedLowerKleenean(LogicalValue::UNLIKELY); } }

    friend Bool same(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return x1._u==x2._u; }
    friend Bool refines(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return x1._u <= x2._u; }
    friend UpperBound<F> refinement(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(min(x1._u,x2._u)); }

    friend Integer cast_integer(UpperBound<F> const& x) {
        return ceil(static_cast<Dyadic>(x._u)); }

    friend OutputStream& operator<<(OutputStream& os, UpperBound<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Bounds<F>::output_places},upward); }
    friend InputStream& operator>>(InputStream& is, UpperBound<F>& x) {
        ARIADNE_NOT_IMPLEMENTED; }
  public:
    friend LowerBound<F> neg(UpperBound<F> const& x);
    friend LowerBound<F> sub(LowerBound<F> const& x1, UpperBound<F> const& x2);
  public:
    friend UpperBound<F> operator*(PositiveBounds<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(mul(up,x2.raw()>=0?x1.upper().raw():x1.lower().raw(),x2.raw())); }
    friend UpperBound<F> operator*(UpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return UpperBound<F>(mul(up,x1.raw(),x1.raw()>=0?x2.upper().raw():x2.lower().raw())); }
    friend UpperBound<F> operator/(UpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return UpperBound<F>(div(up,x1.raw(),x1.raw()>=0?x2.lower().raw():x2.upper().raw())); }
    // Needed to prevent ambiguity; useful as implementation is easier than PositiveBounds version.
    friend UpperBound<F> operator*(PositiveValue<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend UpperBound<F> operator*(UpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return UpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend UpperBound<F> operator/(UpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return UpperBound<F>(div(up,x1.raw(),x2.raw())); }
  private: public:
    static Nat output_places;
    RawType _u;
};

template<class F> template<class FE> UpperBound<F>::UpperBound(Ball<F,FE> const& x) : UpperBound<F>(x.upper_raw()) { }

template<class PR> UpperBound(ValidatedUpperNumber, PR) -> UpperBound<RawFloatType<PR>>;
template<class F> UpperBound(F) -> UpperBound<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(UpperBound<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatUpperBound<PR> FloatFactory<PR>::create(ValidatedUpperNumber const& y) { return FloatUpperBound<PR>(y,_pr); }
template<class PR> inline PositiveFloatUpperBound<PR> FloatFactory<PR>::create(PositiveValidatedUpperNumber const& y) { return PositiveFloatUpperBound<PR>(y,_pr); }

template<class F> class Positive<UpperBound<F>> : public UpperBound<F>
    , DefineConcreteGenericOperators<PositiveUpperBound<F>>
{
    using typename UpperBound<F>::PR;
  public:
    Positive<UpperBound<F>>() : UpperBound<F>() { }
    explicit Positive<UpperBound<F>>(PR const& pr) : UpperBound<F>(pr) { }
    explicit Positive<UpperBound<F>>(F const& x) : UpperBound<F>(x) { }
    template<BuiltinUnsignedIntegral M> Positive<UpperBound<F>>(M m, PR pr) : UpperBound<F>(m,pr) { }
    template<BuiltinUnsignedIntegral M> PositiveValue<F> create(M m) const { return PositiveValue<F>(m,this->precision()); }
    explicit Positive<UpperBound<F>>(UpperBound<F> const& x) : UpperBound<F>(x) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"x="<<x); }
    Positive<UpperBound<F>>(PositiveValidatedUpperNumber const& y, PR pr) : UpperBound<F>(y,pr) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"y="<<y); }
    Positive<UpperBound<F>>(PositiveValue<F> const& x) : UpperBound<F>(x) { }
    Positive<UpperBound<F>>(PositiveBounds<F> const& x) : UpperBound<F>(x) { }
  public:
    friend PositiveUpperBound<F> nul(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(nul(x.raw())); }
    friend PositiveUpperBound<F> pos(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(pos(x.raw())); }
    friend PositiveUpperBound<F> sqr(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(sqr(up,x.raw())); }
    friend PositiveUpperBound<F> rec(PositiveLowerBound<F> const& x) {
        return PositiveUpperBound<F>(rec(up,x.raw())); }
    friend PositiveUpperBound<F> add(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(add(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> mul(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> div(PositiveUpperBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> max(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(max(x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> max(PositiveUpperBound<F> const& x1, UpperBound<F> const& x2) {
        return PositiveUpperBound<F>(max(x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> max(UpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(max(x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> min(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(min(x1.raw(),x2.raw())); }

    friend PositiveLowerBound<F> rec(PositiveUpperBound<F> const& x);
    friend PositiveLowerBound<F> div(PositiveLowerBound<F> const& x1, PositiveUpperBound<F> const& x2);

    friend PositiveUpperBound<F> pow(PositiveUpperBound<F> const& x, Nat m) {
        return PositiveUpperBound<F>(pow(up,x._u,static_cast<Int>(m))); }
    friend Approximation<F> pow(UpperBound<F> const& x, Int n) {
        return pow(Approximation<F>(x),n); }
    // FIXME: Implement pow for Natural/Integer
    friend PositiveUpperBound<F> pow(PositiveUpperBound<F> const& x, Natural const& m) {
        return PositiveUpperBound<F>(pow(up,x._u,static_cast<Int>(m.get_si()))); }

    friend PositiveUpperBound<F> sqrt(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(sqrt(up,x.raw())); }
    friend PositiveUpperBound<F> exp(UpperBound<F> const& x) {
        return PositiveUpperBound<F>(exp(up,x.raw())); }
    friend UpperBound<F> log(PositiveUpperBound<F> const& x) {
        return UpperBound<F>(log(up,x.raw())); }
    friend UpperBound<F> atan(UpperBound<F> const& x) {
        return UpperBound<F>(atan(up,x.raw())); }
    friend PositiveUpperBound<F> atan(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(atan(up,x.raw())); }

  public:
    friend PositiveUpperBound<F> operator+(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(add(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator/(PositiveUpperBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator*=(PositiveUpperBound<F>& x1, PositiveUpperBound<F> const& x2) {
        return x1=x1*x2; }
  public:
    // Needed to prevent ambiguity
    friend PositiveUpperBound<F> operator*(PositiveBounds<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.upper().raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.upper().raw())); }
    friend PositiveUpperBound<F> operator/(PositiveUpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),x2.lower().raw())); }
    friend PositiveUpperBound<F> operator*(PositiveValue<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator/(PositiveUpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator/(PositiveValue<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveUpperBound<F>(div(down,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator/(PositiveUpperBound<F> const& x1, Nat m2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),F(m2,x1.precision()))); }
    friend PositiveLowerBound<F> operator/(PositiveValue<F> const& x1, PositiveUpperBound<F> const& x2);
};

template<class F> inline PositiveUpperBound<F> cast_positive(UpperBound<F> const& x) {
    return PositiveUpperBound<F>(x); }

}

#endif
