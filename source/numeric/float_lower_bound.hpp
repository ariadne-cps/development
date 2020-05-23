/***************************************************************************
 *            numeric/float_lower_bound.hpp
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

/*! \file numeric/float_lower_bound.hpp
 *  \brief Floating-point lower bounds for real numbers.
 */

#ifndef FLOAT_LOWER_BOUND_H
#define FLOAT_LOWER_BOUND_H

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
//! \sa LowerReal, FloatDP, FloatMP, Bounds, UpperBound.
template<class F> class LowerBound
    : public DefineDirectedGroupOperators<LowerBound<F>,UpperBound<F>>
    , public DefineDirectedGroupOperators<UpperBound<F>,LowerBound<F>>
    , public DefineDirectedComparisonOperators<LowerBound<F>,UpperBound<F>,LessTrait<LowerBound<F>>,EqualsTrait<LowerBound<F>>>
    , public DefineDirectedComparisonOperators<UpperBound<F>,LowerBound<F>,LessTrait<UpperBound<F>>,EqualsTrait<UpperBound<F>>>
    , public DefineConcreteGenericOperators<LowerBound<F>>
    , public DeclareFloatOperations<Approximation<F>>
{
  protected:
    typedef LowerTag P; typedef typename F::RoundingModeType RND; typedef typename F::PrecisionType PR;
  public:
    //! <p/>
    typedef LowerTag Paradigm;
    //! <p/>
    typedef LowerBound<F> NumericType;
    //! <p/>
    typedef ValidatedLowerNumber GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR PropertiesType;
  public:
    //! A lower bound of zero with precision \a pr.
    explicit LowerBound<F>(PrecisionType pr) : _l(0.0_x,pr) { }
    //! A lower bound with value \a l.
    explicit LowerBound<F>(RawType const& l) : _l(l) { }

    template<BuiltinIntegral N> LowerBound<F>(N n, PR pr) : LowerBound<F>(ExactDouble(n),pr) { }
    LowerBound<F>(const ExactDouble& d, PR pr) : _l(d,pr) { }
        LowerBound<F>(const TwoExp& t, PR pr) : _l(t,pr) { }
        LowerBound<F>(const Integer& z, PR pr) : _l(z,down,pr) { }
        LowerBound<F>(const Dyadic& w, PR pr) : _l(w,down,pr) { }
        LowerBound<F>(const Decimal& d, PR pr) : _l(d,down,pr) { }
        LowerBound<F>(const Rational& q, PR pr) : _l(q,down,pr) { }
        LowerBound<F>(const Real& r, PR pr);
        LowerBound<F>(const Value<F>& x, PR pr); // FIXME: Should not be necessary
        LowerBound<F>(const Bounds<F>& x, PR pr);
    LowerBound<F>(const LowerBound<F>& x, PR pr);
    //! A lower bound of type \p F from a generic lower bound \a y.
    LowerBound<F>(const ValidatedLowerNumber& y, PR pr);
    template<class FF> requires Constructible<F,FF,RND,PR>
        LowerBound<F>(const LowerBound<FF>& x, PR pr) : _l(x.raw(),down,pr) { }

    //! Convert from lower \em and upper bounds on a number.
    LowerBound<F>(Bounds<F> const& x);
    template<class FE> LowerBound<F>(Ball<F,FE> const& x);
    LowerBound<F>(Value<F> const& x);

        LowerBound<F>& operator=(const Value<F>& x) { return *this=LowerBound<F>(x); }
    //! Assign from the lower bound \a y, keeping the same properties.
    LowerBound<F>& operator=(const ValidatedLowerNumber& y) { return *this=LowerBound<F>(y,this->precision()); }
    //! Create a lower bound from the generic lower bound \a y with the same properties as \a this.
    LowerBound<F> create(const ValidatedLowerNumber& y) const { return LowerBound<F>(y,this->precision()); }
    //! Create a upper bound from the generic upper bound \a y with the same properties as \a this.
    UpperBound<F> create(const ValidatedUpperNumber& y) const { return UpperBound<F>(y,this->precision()); }

    //! Downcast to a generic lower bound.
    operator ValidatedLowerNumber () const;

    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _l.precision(); }
    //! The compuational properties needed to create the lower bound; equivalent to the precision.
    PropertiesType properties() const { return _l.precision(); }
    //! Downcast to a generic lower bound.
    GenericType generic() const { return this->operator GenericType(); }
    //! The raw data used to represent the lower bound.
    RawType const& raw() const { return _l; }
    //! A mutable reference to the raw data used to represent the lower bound.
    RawType& raw() { return _l; }
    //! Under-approximate by a builtin double-precision value. DEPRECATED \deprecated
    double get_d() const { return _l.get_d(); }
  public: // To be removed
    friend Bool is_nan(LowerBound<F> const& x) {
        return is_nan(x._l); }

    friend LowerBound<F> max(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(max(x1._l,x2._l)); }
    friend LowerBound<F> min(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(min(x1._l,x2._l)); }
    friend Approximation<F> abs(LowerBound<F> const& x) {
        return abs(Approximation<F>(x)); }

    friend LowerBound<F> nul(LowerBound<F> const& x) {
        return LowerBound<F>(pos(x._l)); }
    friend LowerBound<F> pos(LowerBound<F> const& x) {
        return LowerBound<F>(pos(x._l)); }
    friend UpperBound<F> neg(LowerBound<F> const& x) {
        return UpperBound<F>(neg(x._l)); }
    friend LowerBound<F> hlf(LowerBound<F> const& x) {
        return LowerBound<F>(hlf(x._l)); }

    friend LowerBound<F> add(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(add(down,x1._l,x2._l)); }
//    friend Approximation<F> sub(LowerBound<F> const& x1, LowerBound<F> const& x2) {
//        return UpperBound<F>(sub(near,x1._l,x2._l)); }
    friend LowerBound<F> sub(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        return LowerBound<F>(sub(down,x1._l,x2._u)); }

    friend LowerBound<F> sqrt(LowerBound<F> const& x) {
        return LowerBound<F>(sqrt(down,x.raw())); }
    friend LowerBound<F> exp(LowerBound<F> const& x) {
        return LowerBound<F>(exp(down,x.raw())); }
    friend LowerBound<F> log(LowerBound<F> const& x) {
        return LowerBound<F>(log(down,x.raw())); }
    friend LowerBound<F> atan(LowerBound<F> const& x) {
        return LowerBound<F>(atan(down,x.raw())); }

    friend ValidatedNegatedSierpinskian eq(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        if(x1._l>x2._u) { return false; }
        else { return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE); } }
    friend ValidatedUpperKleenean lt(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        if(x1._l>=x2._u) { return false; }
        else { return ValidatedUpperKleenean(LogicalValue::LIKELY); } }

    friend Bool same(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return x1._l==x2._l; }
    friend Bool refines(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return x1._l>=x2._l; }
    friend LowerBound<F> refinement(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(max(x1._l,x2._l)); }

    friend Integer cast_integer(LowerBound<F> const& x) {
        return floor(static_cast<Dyadic>(x._l)); }

    friend OutputStream& operator<<(OutputStream& os, LowerBound<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Bounds<F>::output_places},downward); }
    friend InputStream& operator>>(InputStream& is, LowerBound<F>& x) {
        ARIADNE_NOT_IMPLEMENTED; }
  public:
    friend UpperBound<F> neg(LowerBound<F> const& x);
    friend UpperBound<F> sub(UpperBound<F> const& x1, LowerBound<F> const& x2);
  public:
    friend LowerBound<F> operator*(PositiveBounds<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(mul(down,x2.raw()>=0?x1.lower().raw():x1.upper().raw(),x2.raw())); }
    friend LowerBound<F> operator*(LowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return LowerBound<F>(mul(down,x1.raw(),x1.raw()>=0?x2.lower().raw():x2.upper().raw())); }
    friend LowerBound<F> operator/(LowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return LowerBound<F>(div(down,x1.raw(),x1.raw()>=0?x2.upper().raw():x2.lower().raw())); }
    // Needed to prevent ambiguity; useful as implementation is easier than PositiveBounds version.
    friend LowerBound<F> operator*(PositiveValue<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend LowerBound<F> operator*(LowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return LowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend LowerBound<F> operator/(LowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return LowerBound<F>(div(down,x1.raw(),x2.raw())); }
  private: public:
    static Nat output_places;
    RawType _l;
};

template<class F> template<class FE> LowerBound<F>::LowerBound(Ball<F,FE> const& x) : LowerBound<F>(x.lower_raw()) { }

template<class PR> LowerBound(ValidatedLowerNumber, PR) -> LowerBound<RawFloatType<PR>>;
template<class F> LowerBound(F) -> LowerBound<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(LowerBound<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatLowerBound<PR> FloatFactory<PR>::create(ValidatedLowerNumber const& y) { return FloatLowerBound<PR>(y,_pr); }
template<class PR> inline PositiveFloatLowerBound<PR> FloatFactory<PR>::create(PositiveValidatedLowerNumber const& y) { return PositiveFloatLowerBound<PR>(y,_pr); }

template<class F> class Positive<LowerBound<F>> : public LowerBound<F>
    , DefineConcreteGenericOperators<PositiveLowerBound<F>>
{
    using typename LowerBound<F>::PR;
  public:
    Positive<LowerBound<F>>() : LowerBound<F>() { }
    template<BuiltinUnsignedIntegral M>
        Positive<LowerBound<F>>(M m, PR const& pr) : LowerBound<F>(m,pr) { }
    explicit Positive<LowerBound<F>>(PR const& pr) : LowerBound<F>(pr) { }
    explicit Positive<LowerBound<F>>(F const& x) : LowerBound<F>(x) { }
    explicit Positive<LowerBound<F>>(LowerBound<F> const& x) : LowerBound<F>(x) { }
    Positive<LowerBound<F>>(PositiveValidatedLowerNumber const& y, PR pr) : LowerBound<F>(y,pr) { }
    Positive<LowerBound<F>>(PositiveValue<F> const& x) : LowerBound<F>(x) { }
    Positive<LowerBound<F>>(PositiveBounds<F> const& x) : LowerBound<F>(x) { }
  public:
    friend PositiveLowerBound<F> nul(PositiveLowerBound<F> const& x) {
        return PositiveLowerBound<F>(nul(x.raw())); }
    friend PositiveLowerBound<F> pos(PositiveLowerBound<F> const& x) {
        return PositiveLowerBound<F>(pos(x.raw())); }
    friend PositiveLowerBound<F> sqr(PositiveLowerBound<F> const& x) {
        return PositiveLowerBound<F>(sqr(down,x.raw())); }
    friend PositiveLowerBound<F> rec(PositiveUpperBound<F> const& x) {
        return PositiveLowerBound<F>(rec(down,x.raw())); }
    friend PositiveLowerBound<F> add(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(add(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> mul(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> div(PositiveLowerBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> max(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(max(x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> max(PositiveLowerBound<F> const& x1, LowerBound<F> const& x2) {
        return PositiveLowerBound<F>(max(x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> max(LowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(max(x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> min(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(min(x1.raw(),x2.raw())); }

    friend PositiveUpperBound<F> rec(PositiveLowerBound<F> const& x);
    friend PositiveUpperBound<F> div(PositiveUpperBound<F> const& x1, PositiveLowerBound<F> const& x2);

    friend PositiveLowerBound<F> pow(PositiveLowerBound<F> const& x, Nat m) {
        return PositiveLowerBound<F>(pow(down,x.raw(),static_cast<int>(m))); }
    friend Approximation<F> pow(LowerBound<F> const& x, Int n) {
        return pow(Approximation<F>(x),n); }
  public:
    friend PositiveLowerBound<F> operator+(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(add(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator*(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator/(PositiveLowerBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator*=(PositiveLowerBound<F>& x1, PositiveLowerBound<F> const& x2) {
        return x1=x1*x2; }
  public:
    // Needed to prevent ambiguity
    friend PositiveLowerBound<F> operator*(PositiveBounds<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.lower().raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator*(PositiveLowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.lower().raw())); }
    friend PositiveLowerBound<F> operator/(PositiveLowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),x2.upper().raw())); }
    friend PositiveLowerBound<F> operator*(PositiveValue<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator*(PositiveLowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator/(PositiveLowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator/(PositiveValue<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator/(PositiveLowerBound<F> const& x1, Nat m2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),F(m2,x1.precision()))); }
    friend PositiveUpperBound<F> operator/(PositiveValue<F> const& x1, PositiveLowerBound<F> const& x2);
};

template<class F> inline PositiveLowerBound<F> cast_positive(LowerBound<F> const& x) {
    return PositiveLowerBound<F>(x); }

}

#endif
