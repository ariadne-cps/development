/***************************************************************************
 *            numeric/float_value.hpp
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

/*! \file numeric/float_value.hpp
 *  \brief Exact Floating-point representations of real numbers.
 */

#ifndef ARIADNE_FLOAT_VALUE_HPP
#define ARIADNE_FLOAT_VALUE_HPP

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"

#include "integer.hpp"
#include "dyadic.hpp"

namespace Ariadne {

struct DefaultTag;

static_assert(not IsGenericNumericType<FloatValue<DoublePrecision>>::value,"");
static_assert(not IsGenericNumericType<FloatValue<MultiplePrecision>>::value,"");

extern const FloatDPValue infty;

FloatDPValue operator"" _exact(long double lx);

//! \ingroup NumericModule
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
//! \sa Dyadic, Real, FloatDP, FloatMP, Ball, Bounds, Approximation.
template<class F> class Value
    : public DefineFieldOperators<Value<F>,Bounds<F>>
    , public DefineComparisonOperators<Value<F>,LessTrait<Value<F>>,EqualsTrait<Value<F>>>
    , public DefineConcreteGenericOperators<Value<F>>
    , DefineMixedComparisonOperators<Value<F>,ExactNumber,Boolean>
    , DefineMixedComparisonOperators<Value<F>,Rational,Boolean>
    , DefineMixedComparisonOperators<Value<F>,Dyadic,Boolean>
{
  protected:
    typedef ExactTag P; typedef typename F::PrecisionType PR;
  public:
    //! <p/>
    typedef ExactTag Paradigm;
    //! <p/>
    typedef Value<F> NumericType;
    //! <p/>
    typedef ExactNumber GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR PropertiesType;
  public:
    //! Construct an exact numeric value with precision \a pr.
    explicit Value<F>(PrecisionType pr) : _v(0.0_x,pr) { }
    //! Treat \a v as an exact numeric value.
    explicit Value<F>(RawType const& v) : _v(v) { }

    template<BuiltinIntegral N> Value<F>(N n, PR pr) : Value<F>(Integer(n),pr) { }
    Value<F>(const ExactDouble& d, PR pr);
    Value<F>(const TwoExp& t, PR pr);
    Value<F>(const Integer& z, PR pr);
    //! Construct a floating-point value with precision \a pr exactly from the dyadic number \a w.
    //! \pre Requires that \a w can be exactly represented by a value of type \p F with precision \a pr.
    Value<F>(const Dyadic& w, PR pr);
    Value<F>(const Value<F>& x, PR pr);

    template<BuiltinIntegral N> Value<F>& operator=(N n) { _v=n; return *this; }
    Value<F>& operator=(const Integer& z);
    Value<F>& operator=(const TwoExp& t);
    //! Assign from a dyadic number \a w, keeping the same precision.
    //! \pre Requires that \a w can be exactly represented by a value of type \p F with the current precision.
    Value<F>& operator=(const Dyadic& w);

    //! Convert to a generic exact number.
    operator ExactNumber () const;
    //! Convert to a dyadic number.
    explicit operator Dyadic () const;
    explicit operator Rational () const;

    Ball<F> create(ValidatedNumber const&) const;

    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _v.precision(); }
    //! The compuational properties needed to create the bounds; equivalent to the precision.
    PropertiesType properties() const { return _v.precision(); }
    //! Downcast to generic validated bounds.
    GenericType generic() const { return this->operator GenericType(); }
    //! Create a value with the value at its centre, and error bound \a e.
    template<class FE> Ball<F,FE> pm(Error<FE> const& e) const;
    //! The raw data used to represent the value.
    RawType const& raw() const { return _v; }
    //! A mutable reference to the raw data used to represent the value.
    RawType& raw() { return _v; }
    //! Approximate by a builtin double-precision value. DEPRECATED
    double get_d() const { return _v.get_d(); }
  public:
    friend Bool is_nan(Value<F> const& x) {
        return is_nan(x._v); }

    friend Value<F> max(Value<F> const& x1,  Value<F> const& x2) {
        return Value<F>(max(x1._v,x2._v)); }
    friend Value<F> min(Value<F> const& x1,  Value<F> const& x2) {
        return Value<F>(min(x1._v,x2._v)); }
    friend Value<F> abs(Value<F> const& x) {
        return Value<F>(abs(x._v)); }
    friend PositiveLowerBound<F> mig(Value<F> const& x) {
        return PositiveLowerBound<F>(abs(x._v)); }
    friend PositiveUpperBound<F> mag(Value<F> const& x) {
        return PositiveUpperBound<F>(abs(x._v)); }

    friend Value<F> nul(Value<F> const& x) {
        return Value<F>(nul(x._v)); }
    friend Value<F> pos(Value<F> const& x) {
        return Value<F>(pos(x._v)); }
    friend Value<F> neg(Value<F> const& x) {
        return Value<F>(neg(x._v)); }
    friend Value<F> hlf(Value<F> const& x) {
        return Value<F>(hlf(x._v)); }
    friend Bounds<F> sqr(Value<F> const& x) {
        return Bounds<F>(mul(down,x._v,x._v),mul(up,x._v,x._v)); }
    friend Bounds<F> rec(Value<F> const& x) {
        return Bounds<F>(rec(down,x._v),rec(up,x._v)); }
    friend Bounds<F> add(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(add(down,x1._v,x2._v),add(up,x1._v,x2._v)); }
    friend Bounds<F> sub(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(sub(down,x1._v,x2._v),sub(up,x1._v,x2._v)); }
    friend Bounds<F> mul(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(mul(down,x1._v,x2._v),mul(up,x1._v,x2._v)); }
    friend Bounds<F> div(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(div(down,x1._v,x2._v),div(up,x1._v,x2._v)); }

    template<class PRE> friend Ball<F,RawFloatType<PRE>> add(Value<F> const& x1, Value<F> const& x2, PRE pre) {
        typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_add(x1,x2,pre); }
    template<class PRE> friend Ball<F,RawFloatType<PRE>> sub(Value<F> const& x1, Value<F> const& x2, PRE pre) {
        typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_sub(x1,x2,pre); }
    template<class PRE> friend Ball<F,RawFloatType<PRE>> mul(Value<F> const& x1, Value<F> const& x2, PRE pre) {
        typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_mul(x1,x2,pre); }
    template<class PRE> friend Ball<F,RawFloatType<PRE>> div(Value<F> const& x1, Value<F> const& x2, PRE pre) {
        typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_div(x1,x2,pre); }

    template<class FE> friend Value<F> add(Value<F> const& x1, Value<F> const& x2, Error<FE>& e);
    template<class FE> friend Value<F> sub(Value<F> const& x1, Value<F> const& x2, Error<FE>& e);
    template<class FE> friend Value<F> mul(Value<F> const& x1, Value<F> const& x2, Error<FE>& e);
    template<class FE> friend Value<F> div(Value<F> const& x1, Value<F> const& x2, Error<FE>& e);

    friend Value<F> mul(Value<F> const& x, TwoExp const& y) {
        Value<F> yv(y,x.precision()); return Value<F>(mul(near,x.raw(),yv.raw())); }
    friend Value<F> div(Value<F> const& x, TwoExp const& y) {
        Value<F> yv(y,x.precision()); return Value<F>(div(near,x.raw(),yv.raw())); }

    friend Bounds<F> pow(Value<F> const& x, Nat m) {
        return pow(Bounds<F>(x),m); }
    friend Bounds<F> pow(Value<F> const& x, Int n) {
        return pow(Bounds<F>(x),n); }

    friend Bounds<F> med(Value<F> const& x1, Value<F> const& x2) {
        return add(hlf(x1),hlf(x2)); }
    friend Bounds<F> rad(Value<F> const& x1, Value<F> const& x2) {
        return sub(hlf(x2),hlf(x1)); }

    friend Bounds<F> sqrt(Value<F> const& x) {
        return Bounds<F>(sqrt(down,x._v),sqrt(up,x._v)); }
    friend Bounds<F> exp(Value<F> const& x) {
        return Bounds<F>(exp(down,x._v),exp(up,x._v)); }
    friend Bounds<F> log(Value<F> const& x) {
        return Bounds<F>(log(down,x._v),log(up,x._v)); }
    friend Bounds<F> sin(Value<F> const& x) {
        return sin(Bounds<F>(x)); }
    friend Bounds<F> cos(Value<F> const& x) {
        return cos(Bounds<F>(x)); }
    friend Bounds<F> tan(Value<F> const& x) {
        return tan(Bounds<F>(x)); }
    friend Bounds<F> asin(Value<F> const& x) {
        return asin(Bounds<F>(x)); }
    friend Bounds<F> acos(Value<F> const& x) {
        return acos(Bounds<F>(x)); }
    friend Bounds<F> atan(Value<F> const& x) {
        return Bounds<F>(atan(down,x._v),atan(up,x._v)); }

    friend Boolean eq(Value<F> const& x1, Value<F> const& x2) {
        return x1._v == x2._v; }
    friend Boolean lt(Value<F> const& x1, Value<F> const& x2) {
        return x1._v <  x2._v; }

    //! <p/>
    friend Bool same(Value<F> const& x1, Value<F> const& x2) {
        return x1._v==x2._v; }

    //! <p/>
    friend OutputStream& operator<<(OutputStream& os, Value<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Value<F>::output_places},to_nearest); }
    //! <p/>
    friend InputStream& operator>>(InputStream& is, Value<F>& x) {
        auto v = nul(x._v); is >> v; ARIADNE_ASSERT(not is.fail()); x._v=v; return is;}

    friend Integer cast_integer(Value<F> const& x) {
        Dyadic w(x); Integer z=round(w); ARIADNE_ASSERT_MSG(z==w,"Cannot cast non-integral value "<<z<<" to an Integer"); return z; }
  public:
    friend Value<F> max(Value<F> const& x1, Dyadic const& w2) {
        return max(x1,Value<F>(w2,x1.precision())); }
    friend Value<F> max(Dyadic const& w1, Value<F> const& x2) {
        return max(Value<F>(w1,x2.precision()),x2); }
    friend Value<F> min(Value<F> const& x1, Dyadic const& w2) {
        return min(x1,Value<F>(w2,x1.precision())); }
    friend Value<F> min(Dyadic const& w1, Value<F> const& x2) {
        return min(Value<F>(w1,x2.precision()),x2); }
  public:
    friend Value<F> operator*(TwoExp const&, Value<F> const&);
    friend Value<F> operator*(Value<F> const&, TwoExp const&);
    friend Value<F> operator/(Value<F> const&, TwoExp const&);
    friend Value<F>& operator*=(Value<F>&, TwoExp const&);
    friend Value<F>& operator/=(Value<F>&, TwoExp const&);
  public:
    friend Comparison cmp(Value<F> const& x1, Rational const& q2) { return cmp(x1.raw(),q2); }
    friend Comparison cmp(Value<F> const& x1, Dyadic const& w2) { return cmp(x1.raw(),w2); }
    friend Comparison cmp(Value<F> const& x1, Integer const& z2) { return cmp(x1.raw(),z2); }
  public:
    static Nat output_places;
    //! Set the number of decimal places used for the output. DEPRECATED
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawType _v;
  private:
    friend Value<F> shft(Value<F> const& x, Int n) {
        return Value<F>(shft(x.raw(),n)); }
    friend Value<F> operator*(TwoExp const& y, Value<F> const& x) {
        return Value<F>(mul(near,F(y,x.precision()),x.raw())); }
    friend Value<F> operator*(Value<F> const& x, TwoExp const& y) {
        return Value<F>(mul(near,x.raw(),F(y,x.precision()))); }
    friend Value<F> operator/(Value<F> const& x, TwoExp const& y) {
        return Value<F>(div(near,x.raw(),F(y,x.precision()))); }
    friend Value<F>& operator*=(Value<F>& x, TwoExp const& y) { return x=x*y; }
    friend Value<F>& operator/=(Value<F>& x, TwoExp const& y) { return x=x/y; }
};

template<class PR> Value(Dyadic, PR) -> Value<RawFloatType<PR>>;
template<class F> Value(F) -> Value<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(Value<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Dyadic const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Integer const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(ExactDouble const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<BuiltinSignedIntegral N> inline
FloatValue<PR> FloatFactory<PR>::create(N const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<BuiltinUnsignedIntegral M> inline
PositiveFloatValue<PR> FloatFactory<PR>::create(M const& y) { return PositiveFloatValue<PR>(y,_pr); }

template<class F> class Positive<Value<F>> : public Value<F> {
    using typename Value<F>::PR;
  public:
    Positive<Value<F>>() : Value<F>() { }
    explicit Positive<Value<F>>(PR const& pr) : Value<F>(pr) { }
    template<BuiltinUnsignedIntegral M> Positive<Value<F>>(M m, PR pr) : Value<F>(m,pr) { }
    Positive<Value<F>>(TwoExp const& ex, PR pr) : Value<F>(ex,pr) { }
    explicit Positive<Value<F>>(Dyadic const& w, PR pr) : Value<F>(w,pr) { }
    explicit Positive<Value<F>>(F const& x) : Value<F>(x) { }
    explicit Positive<Value<F>>(Value<F> const& x) : Value<F>(x) { }
  public:
    friend PositiveBounds<F> operator+(PositiveValue<F> const& v1, PositiveValue<F> const& v2) {
        return cast_positive(static_cast<Value<F>const&>(v1)+static_cast<Value<F>const&>(v2)); }
    friend PositiveBounds<F> operator*(PositiveValue<F> const& v1, PositiveValue<F> const& v2) {
        return cast_positive(static_cast<Value<F>const&>(v1)*static_cast<Value<F>const&>(v2)); }

    friend Positive<Value<F>> nul(Positive<Value<F>> const& x) { return PositiveValue<F>(nul(x._v)); }
    friend Positive<Value<F>> pos(Positive<Value<F>> const& x) { return PositiveValue<F>(pos(x._v)); }
    friend Positive<Value<F>> hlf(Positive<Value<F>> const& x) { return PositiveValue<F>(hlf(x._v)); }
    friend Positive<Bounds<F>> pow(Positive<Value<F>> const& x, Nat m) {
        return pow(Positive<Bounds<F>>(x),m); }
    friend Positive<Bounds<F>> pow(Positive<Value<F>> const& x, Int n) {
        return pow(Positive<Bounds<F>>(x),n); }
};

template<class F> inline PositiveValue<F> cast_positive(Value<F> const& x) {
    return PositiveValue<F>(x); }

static_assert(IsSame<decltype(declval<FloatDPValue>() < declval<Rational>()),Boolean>::value,"");

extern template Ariadne::Nat Ariadne::Value<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Value<Ariadne::FloatMP>::output_places;


}

#endif
