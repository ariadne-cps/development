/***************************************************************************
 *            float_lower_bound.hpp
 *
 *  Copyright 2008-17  LOWER_BOUNDieter Collins
 *
 ****************************************************************************/

/*
 *  LOWER_BOUNDhis program is free software; you can redistribute it and/or modify
 *  it under the terms of the GLOWER_BOUNDU General Public License as published by
 *  the Free LOWER_BOUNDoftware Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  LOWER_BOUNDhis program is distributed in the hope that it will be useful,
 *  but WLOWER_BOUNDTHOUT ANY WARRANTY; without even the implied warranty of
 *  LOWER_BOUNDERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GLOWER_BOUNDU Library General Public License for more details.
 *
 *  You should have received a copy of the GLOWER_BOUNDU General Public License
 *  along with this program; if not, write to the Free LOWER_BOUNDoftware
 *  Foundation, LOWER_BOUNDnc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file float_lower_bound.hpp
 *  \brief Floating-point lower bounds for real numbers.
 */

#ifndef FLOAT_LOWER_BOUND_H
#define FLOAT_LOWER_BOUND_H

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatLowerBound<PR>> {
    typedef ValidatedLowerNumber GenericType;
    typedef FloatUpperBound<PR> OppositeType;
    typedef PositiveFloatLowerBound<PR> PositiveType;
    typedef ValidatedLowerKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
//! \sa UpperReal, Float64, FloatMP, FloatBounds, FloatUpperBound.
template<class PR> class FloatLowerBound
    : public DispatchDirectedFloatOperations<FloatLowerBound<PR>>
    , public DispatchFloatOperations<FloatApproximation<PR>>
{
    typedef LowerTag P; typedef RawFloat<PR> FLT;
  public:
    typedef LowerTag Paradigm;
    typedef FloatLowerBound<PR> NumericType;
    typedef ValidatedLowerNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    FloatLowerBound<PR>() : _l(0.0) { }
    explicit FloatLowerBound<PR>(PrecisionType pr) : _l(0.0,pr) { }
    explicit FloatLowerBound<PR>(RawFloatType const& l) : _l(l) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> FloatLowerBound<PR>(N n, PR pr) : FloatLowerBound<PR>(ExactDouble(n),pr) { }
    FloatLowerBound<PR>(ExactDouble d, PR pr);
        FloatLowerBound<PR>(const Integer& z, PR pr);
        FloatLowerBound<PR>(const Dyadic& w, PR pr);
        FloatLowerBound<PR>(const Decimal& d, PR pr);
        FloatLowerBound<PR>(const Rational& q, PR pr);
        FloatLowerBound<PR>(const Real& r, PR pr);
    FloatLowerBound<PR>(const FloatLowerBound<PR>& x, PR pr);
    FloatLowerBound<PR>(const ValidatedLowerNumber& y, PR pr);

    FloatLowerBound<PR>(FloatBounds<PR> const& x);
    FloatLowerBound<PR>(FloatBall<PR> const& x);
    FloatLowerBound<PR>(FloatValue<PR> const& x);

        FloatLowerBound<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatLowerBound<PR>(x); }
    FloatLowerBound<PR>& operator=(const ValidatedLowerNumber&);
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y) const;
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y) const;

    operator ValidatedLowerNumber () const;

    PrecisionType precision() const { return _l.precision(); }
    PropertiesType properties() const { return _l.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawFloatType const& raw() const { return _l; }
    RawFloatType& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  public: // To be removed
    friend Bool same(FloatLowerBound<PR> const&, FloatLowerBound<PR> const&);
    friend Bool refines(FloatLowerBound<PR> const&, FloatLowerBound<PR> const&);
    friend FloatLowerBound<PR> refinement(FloatLowerBound<PR> const&, FloatLowerBound<PR> const&);
  public:
    friend FloatLowerBound<PR> operator*(FloatLowerBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return FloatLowerBound<PR>(mul(down,x1.raw(),x2.raw())); }
    friend FloatLowerBound<PR> operator*(FloatLowerBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return FloatLowerBound<PR>(mul(down,x1.raw(),x1.raw()>=0?x2.lower().raw():x2.upper().raw())); }
    friend FloatLowerBound<PR> operator/(FloatLowerBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return FloatLowerBound<PR>(div(down,x1.raw(),x2.raw())); }
    friend FloatLowerBound<PR> operator/(FloatLowerBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return FloatLowerBound<PR>(div(down,x1.raw(),x1.raw()>=0?x2.upper().raw():x2.lower().raw())); }
  private: public:
    static Nat output_places;
    RawFloatType _l;
};

template<class PR> inline FloatFactory<PR> factory(FloatLowerBound<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatLowerBound<PR> FloatFactory<PR>::create(Number<LowerTag> const& y) { return FloatLowerBound<PR>(y,_pr); }

template<class PR> class Positive<FloatLowerBound<PR>> : public FloatLowerBound<PR>
    , public DispatchPositiveDirectedFloatOperations<PositiveFloatLowerBound<PR>,PositiveFloatUpperBound<PR>>
{
  public:
    Positive<FloatLowerBound<PR>>() : FloatLowerBound<PR>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<FloatLowerBound<PR>>(M m) : FloatLowerBound<PR>(m) { }
    explicit Positive<FloatLowerBound<PR>>(RawFloat<PR> const& x) : FloatLowerBound<PR>(x) { }
    explicit Positive<FloatLowerBound<PR>>(FloatLowerBound<PR> const& x) : FloatLowerBound<PR>(x) { }
    explicit Positive<FloatLowerBound<PR>>(ValidatedLowerNumber const& y, PR pr) : FloatLowerBound<PR>(y,pr) { }
    Positive<FloatLowerBound<PR>>(PositiveFloatValue<PR> const& x) : FloatLowerBound<PR>(x) { }
    Positive<FloatLowerBound<PR>>(PositiveFloatBounds<PR> const& x) : FloatLowerBound<PR>(x) { }
  public:
    friend PositiveFloatLowerBound<PR> operator*(PositiveFloatLowerBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return PositiveFloatLowerBound<PR>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveFloatLowerBound<PR> operator*(PositiveFloatLowerBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return PositiveFloatLowerBound<PR>(mul(down,x1.raw(),x2.lower().raw())); }
    friend PositiveFloatLowerBound<PR> operator/(PositiveFloatLowerBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return PositiveFloatLowerBound<PR>(div(down,x1.raw(),x2.raw())); }
    friend PositiveFloatLowerBound<PR> operator/(PositiveFloatLowerBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return PositiveFloatLowerBound<PR>(div(down,x1.raw(),x2.upper().raw())); }
};

template<class PR> inline PositiveFloatLowerBound<PR> cast_positive(FloatLowerBound<PR> const& x) {
    return PositiveFloatLowerBound<PR>(x); }

}

#endif