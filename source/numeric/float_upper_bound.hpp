/***************************************************************************
 *            float_upper_bound.hpp
 *
 *  Copyright 2008-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file float_upper_bound.hpp
 *  \brief Floating-point upper bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_UPPER_BOUND_HPP
#define ARIADNE_FLOAT_UPPER_BOUND_HPP

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatUpperBound<PR>> {
    typedef ValidatedUpperNumber GenericType;
    typedef FloatLowerBound<PR> OppositeType;
    typedef PositiveFloatUpperBound<PR> PositiveType;
    typedef ValidatedUpperKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
//! \sa UpperReal, Float64, FloatMP, FloatBounds, FloatLowerBound.
template<class PR> class FloatUpperBound
    : public DispatchDirectedFloatOperations<FloatUpperBound<PR>>
    , public DispatchFloatOperations<FloatApproximation<PR>>
{
    typedef UpperTag P; typedef RawFloat<PR> FLT;
  public:
    typedef UpperTag Paradigm;
    typedef FloatUpperBound<PR> NumericType;
    typedef ValidatedUpperNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    FloatUpperBound<PR>() : _u(0.0) { }
    explicit FloatUpperBound<PR>(PrecisionType pr) : _u(0.0,pr) { }
    explicit FloatUpperBound<PR>(RawFloatType const& u) : _u(u) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> FloatUpperBound<PR>(N n, PR pr) : FloatUpperBound<PR>(ExactDouble(n),pr) { }
    FloatUpperBound<PR>(ExactDouble d, PR pr);
        FloatUpperBound<PR>(const Integer& z, PR pr);
        FloatUpperBound<PR>(const Dyadic& w, PR pr);
        FloatUpperBound<PR>(const Decimal& d, PR pr);
        FloatUpperBound<PR>(const Rational& q, PR pr);
        FloatUpperBound<PR>(const Real& r, PR pr);
    FloatUpperBound<PR>(const FloatUpperBound<PR>& x, PR pr);
    FloatUpperBound<PR>(const ValidatedUpperNumber& y, PR pr);

    FloatUpperBound<PR>(FloatBounds<PR> const& x);
    FloatUpperBound<PR>(FloatBall<PR> const& x);
    FloatUpperBound<PR>(FloatValue<PR> const& x);
    FloatUpperBound<PR>(FloatError<PR> const& x); // FIXME: Remove

        FloatUpperBound<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatUpperBound<PR>(x); }
    FloatUpperBound<PR>& operator=(const ValidatedUpperNumber& y);
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y) const;
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y) const;

    operator ValidatedUpperNumber () const;

    PrecisionType precision() const { return _u.precision(); }
    PropertiesType properties() const { return _u.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawFloatType const& raw() const { return _u; }
    RawFloatType& raw() { return _u; }
    double get_d() const { return _u.get_d(); }
  public: // To be removed
    friend Bool same(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
    friend Bool refines(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
    friend FloatUpperBound<PR> refinement(FloatUpperBound<PR> const&, FloatUpperBound<PR> const&);
  public:
    friend FloatUpperBound<PR> operator*(FloatUpperBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return FloatUpperBound<PR>(mul(up,x1.raw(),x2.raw())); }
    friend FloatUpperBound<PR> operator*(FloatUpperBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return FloatUpperBound<PR>(mul(up,x1.raw(),x1.raw()>=0?x2.upper().raw():x2.lower().raw())); }
    friend FloatUpperBound<PR> operator/(FloatUpperBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return FloatUpperBound<PR>(div(up,x1.raw(),x2.raw())); }
    friend FloatUpperBound<PR> operator/(FloatUpperBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return FloatUpperBound<PR>(div(up,x1.raw(),x1.raw()>=0?x2.lower().raw():x2.upper().raw())); }
  private: public:
    static Nat output_places;
    RawFloatType _u;
};

template<class PR> inline FloatFactory<PR> factory(FloatUpperBound<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatUpperBound<PR> FloatFactory<PR>::create(Number<UpperTag> const& y) { return FloatUpperBound<PR>(y,_pr); }

template<class PR> class Positive<FloatUpperBound<PR>> : public FloatUpperBound<PR>
    , public DispatchPositiveDirectedFloatOperations<PositiveFloatUpperBound<PR>,PositiveFloatLowerBound<PR>>
{
  public:
    Positive<FloatUpperBound<PR>>() : FloatUpperBound<PR>() { }
    explicit Positive<FloatUpperBound<PR>>(PR const& pr) : FloatUpperBound<PR>(pr) { }
    explicit Positive<FloatUpperBound<PR>>(RawFloat<PR> const& x) : FloatUpperBound<PR>(x) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<FloatUpperBound<PR>>(M m, PR pr) : FloatUpperBound<PR>(m,pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> create(M m) const { return PositiveFloatValue<PR>(m,this->precision()); }
    explicit Positive<FloatUpperBound<PR>>(FloatUpperBound<PR> const& x) : FloatUpperBound<PR>(x) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"x="<<x); }
    explicit Positive<FloatUpperBound<PR>>(ValidatedUpperNumber const& y, PR pr) : FloatUpperBound<PR>(y,pr) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"y="<<y); }
    Positive<FloatUpperBound<PR>>(PositiveFloatValue<PR> const& x) : FloatUpperBound<PR>(x) { }
    Positive<FloatUpperBound<PR>>(PositiveFloatBounds<PR> const& x) : FloatUpperBound<PR>(x) { }
  public:
    friend PositiveFloatUpperBound<PR> operator*(PositiveFloatUpperBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return PositiveFloatUpperBound<PR>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveFloatUpperBound<PR> operator*(PositiveFloatUpperBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return PositiveFloatUpperBound<PR>(mul(up,x1.raw(),x2.upper().raw())); }
    friend PositiveFloatUpperBound<PR> operator/(PositiveFloatUpperBound<PR> const& x1, PositiveFloatValue<PR> const& x2) {
        return PositiveFloatUpperBound<PR>(div(up,x1.raw(),x2.raw())); }
    friend PositiveFloatUpperBound<PR> operator/(PositiveFloatUpperBound<PR> const& x1, PositiveFloatBounds<PR> const& x2) {
        return PositiveFloatUpperBound<PR>(div(up,x1.raw(),x2.lower().raw())); }
};

template<class PR> inline PositiveFloatUpperBound<PR> cast_positive(FloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(x); }

}

#endif