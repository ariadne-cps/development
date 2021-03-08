/***************************************************************************
 *            numeric/concepts.hpp
 *
 *  Copyright  2020  Pieter Collins
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

/*! \file numeric/concepts.hpp
 *  \brief
 */

#pragma once

#include "../utility/typedefs.hpp"
#include "../numeric/paradigm.hpp"

#ifndef ARIADNE_NUMERIC_CONCEPTS_HPP
#define ARIADNE_NUMERIC_CONCEPTS_HPP

namespace Ariadne {

enum class Comparison : char;

template<class L> concept IsLogical = Constructible<L,Bool> and requires(L l) {
    { not (not l) } -> ConvertibleTo<L>;
    { l and l } -> ConvertibleTo<L>;
    { l or l } -> ConvertibleTo<L>;
};

template<class X1, class X2=X1> concept HasAdd = requires(X1 x1, X2 x2) { add(x1,x2); };
template<class X1, class X2=X1> concept HasSub = requires(X1 x1, X2 x2) { sub(x1,x2); };
template<class X1, class X2=X1> concept HasMul = requires(X1 x1, X2 x2) { mul(x1,x2); };
template<class X1, class X2=X1> concept HasDiv = requires(X1 x1, X2 x2) { div(x1,x2); };

template<class X> using PosType = decltype(pos(declval<X>()));
template<class X> using NegType = decltype(neg(declval<X>()));
template<class X1, class X2=X1> using AddType = decltype(add(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using SubType = decltype(sub(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using MulType = decltype(mul(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using DivType = decltype(div(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using MaxType = decltype(max(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using MinType = decltype(min(declval<X1>(),declval<X2>()));


template<class Y> concept Lattice = requires(Y y) {
    { abs(y) } -> ConvertibleTo<Y>;
    { max(y,y) } -> SameAs<Y>;
    { min(y,y) } -> SameAs<Y>;
};

template<class Y1, class Y2=Y1> concept Comparible = requires(Y1 y1, Y2 y2) {
    { cmp(y1,y2) } -> SameAs<Comparison>;
};

template<class Y1, class Y2=Y1> concept Ordered = Comparible<Y1,Y2> or requires(Y1 y1, Y2 y2) {
    lt(y1,y2);
    lt(y2,y1);
};

template<class Y1, class Y2=Y1> concept Equality = requires(Y1 y1, Y2 y2) {
    eq(y1,y2);
    eq(y2,y1);
};

template<class Y, class NY> concept DirectedGroups = requires(Y y, NY ny) {
    { neg(y) } -> SameAs<NY>;
    { neg(ny) } -> SameAs<Y>;
    { add(y,y) } -> SameAs<Y>;
    { add(ny,ny) } -> SameAs<NY>;
    { sub(y,ny) } -> SameAs<Y>;
    { sub(ny,y) } -> SameAs<NY>;
};

template<class Y> concept DirectedGroup = DirectedGroups<Y,decltype(neg(declval<Y>()))>;

template<class Y, class QY> concept DirectedSemiFields = requires(Y y, QY qy) {
    { rec(y) } -> SameAs<QY>;
    { rec(qy) } -> SameAs<Y>;

    { add(y,y) } -> SameAs<Y>;
    { add(qy,qy) } -> SameAs<QY>;
    { mul(y,y) } -> SameAs<Y>;
    { mul(qy,qy) } -> SameAs<QY>;
    { div(y,qy) } -> SameAs<Y>;
    { div(qy,y) } -> SameAs<QY>;
};

template<class Y> concept DirectedSemiField = DirectedSemiFields<Y,decltype(rec(declval<Y>()))>;

template<class Y> concept AbelianGroup = (not BuiltinArithmetic<Y>) and requires(Y y) {
    { nul(y) } -> SameAs<Y>;
    { pos(y) } -> SameAs<Y>;
    { neg(y) } -> SameAs<Y>;
//    { add(y1,y2) } -> ConvertibleFrom<Y>;
//    { sub(y1,y2) } -> SameAs<decltype(add(y1,y2))>;
    { add(y,y) } -> SameAs<Y>;
    { sub(y,y) } -> SameAs<Y>;
};

template<class Y> concept Ring = AbelianGroup<Y> and requires(Y y) {
//    { sqr(y,y) } -> ConvertibleTo<AddType<Y,Y>>;
    { mul(y,y) } -> SameAs<AddType<Y,Y>>;
//    { mul(y,y) } -> SameAs<Y>;
//    { fma(y, y, y) } -> SameAs<Y>;
};

template<class Y> concept Field = Ring<Y> and requires(Y y) {
    { hlf(y) } -> SameAs<Y>;
    { rec(y) } -> SameAs<AddType<Y,Y>>;
    { div(y,y) } -> SameAs<AddType<Y,Y>>;
//    { rec(y) } -> SameAs<Y>;
//    { div(y,y) } -> SameAs<Y>;
};

template<class Y> concept OrderedLattice = Ordered<Y> and Lattice<Y>;

template<class Y> concept LatticeRing = Lattice<Y> and Ring<Y>;
template<class Y> concept LatticeField = Lattice<Y> and Field<Y>;

template<class Y> concept OrderedLatticeRing = OrderedLattice<Y> and Ring<Y>;
template<class Y> concept OrderedLatticeField = OrderedLattice<Y> and Field<Y>;


//template<class X, class RND> concept RoundedMode = requires(RND rnd, X x, X x1, X x2) {
template<class X> concept RoundedRing = requires(typename X::RoundingModeType rnd, X x, X x1, X x2) {
    { nul(x) } -> SameAs<X>;
    { pos(x) } -> SameAs<X>;
    { neg(x) } -> SameAs<X>;
    { hlf(x) } -> SameAs<X>;
    { add(rnd, x1, x2) } -> SameAs<X>;
    { sub(rnd, x1, x2) } -> SameAs<X>;
    { mul(rnd, x1, x2) } -> SameAs<X>;
    { fma(rnd, x, x1, x2) } -> SameAs<X>;
};

template<class X> concept RoundedField = RoundedRing<X> and requires(typename X::RoundingModeType rnd, X x, X x1, X x2) {
    { rec(rnd, x) } -> SameAs<X>;
    { div(rnd, x1, x2) } -> SameAs<X>;
};

template<class X> concept RoundedLatticeRing = Lattice<X> and RoundedRing<X>;
template<class X> concept RoundedLatticeField = Lattice<X> and RoundedField<X>;

template<class X> concept IsRoundedRing = RoundedRing<X>;
template<class X> concept IsRoundedField = RoundedField<X>;
template<class X> concept IsRoundedLatticeRing = RoundedLatticeRing<X>;
template<class X> concept IsRoundedLatticeField = RoundedLatticeField<X>;


template<class X1, class X2> concept AreMixedLattice = requires(X1 const& x1, X2 const& x2) {
    max(x1,x2);
    { max(x2,x1) } -> SameAs<decltype(max(x1,x2))>;
    min(x1,x2);
    min(x2,x1);
};

template<class X1, class X2> concept AreMixedRing = requires(X1 const& x1, X2 const& x2) {
    -x1;
    -x2;
    x1+x2;
    { x2+x1 } -> SameAs<decltype(x1+x2)>;
    { x1-x2 } -> SameAs<decltype(x1+(-x2))>;
    { x2-x1 } -> SameAs<decltype(x2+(-x1))>;
    { x1-x2 } -> SameAs<decltype(x1+(-x2))>;
    { x2-x1 } -> SameAs<decltype(x2+(-x1))>;
    x1*x2;
    { x2*x1 } -> SameAs<decltype(x1*x2)>;
};

template<class X1, class X2> concept AreMixedField = AreMixedRing<X1,X2> and requires(X1 const& x1, X2 const& x2) {
    x1/x2;
    x2/x1;
};

template<class X1, class X2> concept AreMixedOrdered = requires(X1 const& x1, X2 const& x2) {
    x1 == x2;
    { x2 == x1 } -> SameAs<decltype(x1==x2)>;
    { x1 != x2 } -> SameAs<decltype(!(x1==x2))>;
    { x2 != x1 } -> SameAs<decltype(x2!=x1)>;
    x1 <= x2;
    x2 <= x1;
    { x1 >= x2 } -> SameAs<decltype(x2<=x1)>;
    { x2 >= x1 } -> SameAs<decltype(x1<=x2)>;
    { x1 <  x2 } -> SameAs<decltype(!(x2<=x1))>;
    { x2 <  x1 } -> SameAs<decltype(!(x1<=x2))>;
    { x1 >  x2 } -> SameAs<decltype(x2< x1)>;
    { x2 >  x1 } -> SameAs<decltype(x2< x1)>;
};

template<class X1, class X2> concept AreMixedOrderedLattice = AreMixedOrdered<X1,X2> and AreMixedLattice<X1,X2>;

template<class X1, class X2> concept AreMixedLatticeRing = AreMixedLattice<X1,X2> and AreMixedRing<X1,X2>;
template<class X1, class X2> concept AreMixedLatticeField = AreMixedLattice<X1,X2> and AreMixedField<X1,X2>;

template<class X1, class X2> concept AreMixedOrderedLatticeRing = AreMixedOrderedLattice<X1,X2> and AreMixedRing<X1,X2>;
template<class X1, class X2> concept AreMixedOrderedLatticeField = AreMixedOrderedLattice<X1,X2> and AreMixedField<X1,X2>;

template<class X> concept IsOrderedLatticeRing = OrderedLatticeRing<X> and AreMixedOrderedLatticeRing<X,X>;

class Real;
class TwoExp;

class DoublePrecision;
class MultiplePrecision;

struct DecimalPlaces;
struct DecimalPrecision;

template<class X> class Positive;

class ApproximateDouble;
class ExactDouble;
class Integer;
class Dyadic;
class Decimal;
class Rational;
class DoublePrecision;
class FloatDP;
struct RoundDownward;
struct RoundUpward;
struct RoundToNearest;
struct RoundApproximately;

template<class P> class Number;
template<class P> class LowerNumber;
template<class P> class UpperNumber;
using ValidatedLowerNumber = LowerNumber<ValidatedTag>;
using ValidatedUpperNumber = UpperNumber<ValidatedTag>;

} // namespace Ariadne

#endif
