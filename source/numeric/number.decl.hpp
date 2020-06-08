/***************************************************************************
 *            numeric/number.decl.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/number.decl.hpp
 *  \brief
 */

#ifndef ARIADNE_NUMBER_DECL_HPP
#define ARIADNE_NUMBER_DECL_HPP

#include <stdexcept>

#include "../utility/metaprogramming.hpp"
#include "../utility/typedefs.hpp"
#include "../numeric/paradigm.hpp"

namespace Ariadne {

class Real;

/************ Number *********************************************************/

class DivideByZeroError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};


template<class X> struct IsNumericType : False { };

template<class X> using NumericType = typename X::NumericType;

template<class X> struct NumericTraits;
template<class X> using GenericTrait = typename NumericTraits<X>::GenericType;
template<class X> using PositiveTrait = typename NumericTraits<X>::PositiveType;
template<class X> using OppositeTrait = typename NumericTraits<X>::OppositeType;
template<class X> using LessTrait = typename NumericTraits<X>::LessType;
template<class X> using EqualsTrait = typename NumericTraits<X>::EqualsType;

template<class X> using PropertiesType = typename X::PropertiesType;
template<class X> using GenericType = typename X::GenericType;

template<class X> using GenericNumericType = GenericType<NumericType<X>>;


template<class T> concept HasGenericType = requires { typename T::GenericType; };
template<class T> concept HasPrecisionType = requires { typename T::PrecisionType; };

template<class X> struct IsConcrete { static const bool value = HasGenericType<X>; };

template<class X> struct IsScalar;
template<class X> struct IsConcreteScalar { static const bool value = IsConcrete<X>::value and IsScalar<X>::value; };
template<class X> struct IsGenericScalar { static const bool value = (not IsConcrete<X>::value) and IsScalar<X>::value; };

template<class X> concept GenericScalar = IsGenericScalar<X>::value;
template<class X> concept ConcreteScalar = IsConcreteScalar<X>::value;

template<class R> struct IsConcreteNumericType : IsConvertible<R,Real> { };
template<class X> concept ConcreteNumber = IsConcreteNumericType<X>::value;


typedef uint Nat;
typedef int Int;
typedef double Dbl;

class Nat32;
class Nat64;
class Int32;
class Int64;

class ApproximateDouble;
class ExactDouble;

class TwoExp;

class Dyadic;
class Decimal;

class Integer;
class Rational;
class Real;

template<class F> class Bounds;
template<> class Bounds<Dyadic>;
using DyadicBounds = Bounds<Dyadic>;

template<class F> class Approximation;
template<> class Approximation<Dyadic>;
using DyadicApproximation = Approximation<Dyadic>;

template<class R, class A> R integer_cast(const A& a);

Integer round(Dyadic const&);
Integer floor(Dyadic const&);
Integer ceil(Dyadic const&);

//! \ingroup NumericModule
//! \brief A modifier declaring that a number is positive.
template<class X> class Positive;

template<class P=Void> class Number;
template<class P> class LowerNumber;
template<class P> class UpperNumber;

template<> struct IsNumericType<Nat>;
template<> struct IsNumericType<Int>;
template<> struct IsNumericType<Dbl>;

template<> struct IsNumericType<Integer>;
template<> struct IsNumericType<Rational>;
template<> struct IsNumericType<Real>;


template<class Y> struct IsGenericNumericType;

template<class Y> struct IsGenericNumericType : IsBuiltinIntegral<Y> { };
template<> struct IsGenericNumericType<ExactDouble> : True { };
template<> struct IsGenericNumericType<Integer> : True { };
template<> struct IsGenericNumericType<Dyadic> : True { };
template<> struct IsGenericNumericType<Rational> : True { };
template<> struct IsGenericNumericType<Real> : True { };
template<> struct IsGenericNumericType<Dbl> : True { };
template<class P> struct IsGenericNumericType<Number<P>> : True { };
template<class P> struct IsGenericNumericType<UpperNumber<P>> : True { };
template<class P> struct IsGenericNumericType<LowerNumber<P>> : True { };
template<class Y> struct IsGenericNumericType<Positive<Y>> : IsGenericNumericType<Y> { };

template<class Y> using IsGenericNumber = IsGenericNumericType<Y>;
template<class Y> concept GenericNumber = IsGenericNumericType<Y>::value;


//! \relates Number
//! \name Type synonyms
//!@{
using ExactNumber=Number<ExactTag>; //!< Alias for generic exact numbers.
using EffectiveNumber=Number<EffectiveTag>; //!< Alias for generic effective numbers.
using EffectiveUpperNumber=UpperNumber<EffectiveTag>; //!< Alias for generic effective upper numbers.
using EffectiveLowerNumber=LowerNumber<EffectiveTag>; //!< Alias for generic effective lower numbers.
using ValidatedNumber=Number<ValidatedTag>; //!< Alias for generic validated numbers.
using ValidatedUpperNumber=UpperNumber<ValidatedTag>; //!< Alias for generic validated upper numbers.
using ValidatedLowerNumber=LowerNumber<ValidatedTag>; //!< Alias for generic validated lower numbers.
using ApproximateNumber=Number<ApproximateTag>; //!< Alias for generic approximate numbers.

template<class P> using PositiveNumber = Positive<Number<P>>; //!< Alias for positive numbers.
template<class P> using PositiveUpperNumber = Positive<UpperNumber<P>>; //!< Alias for positive upper numbers.
template<class P> using PositiveLowerNumber = Positive<LowerNumber<P>>; //!< Alias for positive lower numbers.
using PositiveExactNumber=PositiveNumber<ExactTag>; //!< <p/>
using PositiveEffectiveNumber=PositiveNumber<EffectiveTag>; //!< <p/>
using PositiveEffectiveUpperNumber=PositiveUpperNumber<EffectiveTag>; //!< <p/>
using PositiveEffectiveLowerNumber=PositiveLowerNumber<EffectiveTag>; //!< <p/>
using PositiveValidatedNumber=PositiveNumber<ValidatedTag>; //!< <p/>
using PositiveValidatedUpperNumber=PositiveUpperNumber<ValidatedTag>; //!< <p/>
using PositiveValidatedLowerNumber=PositiveLowerNumber<ValidatedTag>; //!< <p/>
using PositiveApproximateNumber=PositiveNumber<ApproximateTag>; //!< <p/>

using ValidatedErrorNumber = PositiveValidatedUpperNumber; //!< Alias for validated error bounds.
//@!}

} // namespace Ariadne

#endif
