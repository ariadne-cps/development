/***************************************************************************
 *            numeric/float.decl.hpp
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

/*! \file numeric/float.decl.hpp
 *  \brief
 */

#include "../utility/typedefs.hpp"
#include "../numeric/paradigm.hpp"

#ifndef ARIADNE_FLOAT_DECL_HPP
#define ARIADNE_FLOAT_DECL_HPP

namespace Ariadne {

class Real;
class TwoExp;

class DoublePrecision;
class MultiplePrecision;

using DP = DoublePrecision;
using MP = MultiplePrecision;

template<class F> using PrecisionType = typename F::PrecisionType;

struct DecimalPlaces { Nat _places; explicit DecimalPlaces(Nat plc) : _places(plc) { } operator uint() const { return _places; } };
struct DecimalPrecision { Nat _figures; explicit DecimalPrecision(Nat figs) : _figures(figs) { } operator uint() const { return _figures; } };

template<class X> class Positive;

class FloatDP;
class FloatMP;

using DoublePrecisionFloat = FloatDP;
using MultiplePrecisionFloat = FloatMP;

using RawFloatDP = FloatDP;
using RawFloatMP = FloatMP;

#ifdef DOXYGEN
//! \ingroup NumericModule
//! \name Type extraction
//
//!@{
//! \brief The raw floating-point number type corresponding to precision \a PR. \ingroup NumericModule
template<class PR> typedef RawFloatType;
//! \brief The raw floating-point number type corresponding to precision \a PR. \ingroup NumericModule
template<class PR> typedef RawFloat;
//!@}
#else
template<class PR> struct RawFloatTypedef;
template<> struct RawFloatTypedef<DoublePrecision> { typedef FloatDP Type; };
template<> struct RawFloatTypedef<MultiplePrecision> { typedef FloatMP Type; };
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));
template<class PR> using RawFloatType = typename RawFloatTypedef<PR>::Type;

RawFloatDP cast_raw_float(DoublePrecision);
RawFloatMP cast_raw_float(MultiplePrecision);
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));
#endif


template<class F> class Approximation;
template<class F> class LowerBound;
template<class F> class UpperBound;
template<class F> class Bounds;
template<class F, class FE=F> class Ball;
template<class F> class Value;

template<class F> class Error;

//! \relates Positive
//! \name Type shorthands for positive user classes.
//!@{
template<class F> using PositiveApproximation = Positive<Approximation<F>>; //!< <p/>
template<class F> using PositiveLowerBound = Positive<LowerBound<F>>; //!< <p/>
template<class F> using PositiveUpperBound = Positive<UpperBound<F>>; //!< <p/>
template<class F> using PositiveBounds = Positive<Bounds<F>>; //!< <p/>
template<class F, class FE=F> using PositiveBall = Positive<Ball<F>>; //!< <p/>
template<class F> using PositiveValue = Positive<Value<F>>; //!< <p/>
//!@}



//! \ingroup NumericModule
//! \name Type shorthands
//
//!@{
//! \ingroup NumericModule
template<class PR> using FloatApproximation = Approximation<RawFloatType<PR>>; //!< <p/> \ingroup NumericModule
template<class PR> using FloatLowerBound = LowerBound<RawFloatType<PR>>; //!< <p/> \ingroup NumericModule
template<class PR> using FloatUpperBound = UpperBound<RawFloatType<PR>>; //!< <p/> \ingroup NumericModule
template<class PR> using FloatBounds = Bounds<RawFloatType<PR>>; //!< <p/> \ingroup NumericModule
template<class PR, class PRE=PR> using FloatBall = Ball<RawFloatType<PR>,RawFloatType<PRE>>; //!< <p/> \ingroup NumericModule
template<class PR> using FloatValue = Value<RawFloatType<PR>>; //!< <p/> \ingroup NumericModule

template<class PR> using PositiveFloatApproximation = Positive<FloatApproximation<PR>>; //!< <p/> \ingroup NumericModule
template<class PR> using PositiveFloatLowerBound = Positive<FloatLowerBound<PR>>; //!< <p/> \ingroup NumericModule
template<class PR> using PositiveFloatUpperBound = Positive<FloatUpperBound<PR>>; //!< <p/> \ingroup NumericModule
template<class PR> using PositiveFloatBounds = Positive<FloatBounds<PR>>; //!< <p/> \ingroup NumericModule
template<class PR, class PRE=PR> using PositiveFloatBall = Positive<FloatBall<PR,PRE>>; //!< <p/> \ingroup NumericModule
template<class PR> using PositiveFloatValue = Positive<FloatValue<PR>>; //!< <p/> \ingroup NumericModule

template<class PR> using FloatError = Error<RawFloatType<PR>>; //!< <p/> \ingroup NumericModule
//!@}

template<class P, class PR, class PRE=PR> struct FloatTypedef;
template<class PR> struct FloatTypedef<ApproximationTag,PR> { typedef FloatApproximation<PR> Type; };
template<class PR> struct FloatTypedef<LowerTag,PR> { typedef FloatLowerBound<PR> Type; };
template<class PR> struct FloatTypedef<UpperTag,PR> { typedef FloatUpperBound<PR> Type; };
template<class PR> struct FloatTypedef<OrderTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR, class PRE> struct FloatTypedef<MetricTag,PR,PRE> { typedef FloatBall<PR,PRE> Type; };
template<class PR> struct FloatTypedef<ExactTag,PR> { typedef FloatValue<PR> Type; };
template<class PR> struct FloatTypedef<ErrorTag,PR> { typedef FloatError<PR> Type; };

template<class PR> struct FloatTypedef<ApproximateTag,PR> { typedef FloatApproximation<PR> Type; };
template<class PR> struct FloatTypedef<ValidatedTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR> struct FloatTypedef<EffectiveTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR, class PRE> struct FloatTypedef<ValidatedTag,PR,PRE> { typedef FloatBall<PR,PRE> Type; };
template<class PR, class PRE> struct FloatTypedef<EffectiveTag,PR,PRE> { typedef FloatBall<PR,PRE> Type; };

template<class P, class PR, class PRE=PR> using FloatType = typename FloatTypedef<P,PR,PRE>::Type;


//template<class P, class PR, class PRE=PR> using Float = typename FloatTypedef<P,PR,PRE>::Type;
//template<class P> using FloatDP=Float<P,DoublePrecision>;
//template<class P> using FloatMP=Float<P,MultiplePrecision>;

//! \relates FloatDP
//! \name Type synonyms for FloatDP based numbers.
//!@{
using FloatDPApproximation = FloatApproximation<DoublePrecision>; //!< <p/>
using FloatDPLowerBound = FloatLowerBound<DoublePrecision>; //!< <p/>
using FloatDPUpperBound = FloatUpperBound<DoublePrecision>; //!< <p/>
using FloatDPBounds = FloatBounds<DoublePrecision>; //!< <p/>
using FloatDPBall = FloatBall<DoublePrecision>; //!< <p/>
using FloatDPValue = FloatValue<DoublePrecision>; //!< <p/>
using FloatDPError = FloatError<DoublePrecision>; //!< <p/>
using PositiveFloatDPApproximation = PositiveFloatApproximation<DoublePrecision>; //!< <p/>
using PositiveFloatDPLowerBound = PositiveFloatLowerBound<DoublePrecision>; //!< <p/>
using PositiveFloatDPUpperBound = PositiveFloatUpperBound<DoublePrecision>; //!< <p/>
using PositiveFloatDPBounds = PositiveFloatBounds<DoublePrecision>; //!< <p/>
using PositiveFloatDPBall = PositiveFloatBall<DoublePrecision>; //!< <p/>
using PositiveFloatDPValue = PositiveFloatValue<DoublePrecision>; //!< <p/>

using FloatMPDPBall = FloatBall<MultiplePrecision,DoublePrecision>; //!< <p/>
//!@}

void foo();

//! \relates FloatMP
//! \name Type synonyms for FloatMP based numbers.
//!@{
using FloatMPApproximation = FloatApproximation<MultiplePrecision>; //!< \brief <p/>
using FloatMPLowerBound = FloatLowerBound<MultiplePrecision>; //!< <p/>
using FloatMPUpperBound = FloatUpperBound<MultiplePrecision>; //!< <p/>
using FloatMPBounds = FloatBounds<MultiplePrecision>; //!< <p/>
using FloatMPBall = FloatBall<MultiplePrecision>; //!< <p/>
using FloatMPValue = FloatValue<MultiplePrecision>; //!< <p/>
using FloatMPError = FloatError<MultiplePrecision>; //!< <p/>
using PositiveFloatMPApproximation = PositiveFloatApproximation<MultiplePrecision>; //!< <p/>
using PositiveFloatMPLowerBound = PositiveFloatLowerBound<MultiplePrecision>; //!< <p/>
using PositiveFloatMPUpperBound = PositiveFloatUpperBound<MultiplePrecision>; //!< <p/>
using PositiveFloatMPBounds = PositiveFloatBounds<MultiplePrecision>; //!< <p/>
using PositiveFloatMPBall = PositiveFloatBall<MultiplePrecision>; //!< <p/>
using PositiveFloatMPValue = PositiveFloatValue<MultiplePrecision>; //!< <p/>

using FloatMDPBall = FloatBall<MultiplePrecision,DoublePrecision>; //!< <p/>
//!@}


template<class X> struct IsFloat : False { };
template<> struct IsFloat<FloatDP> : True { };
template<> struct IsFloat<FloatMP> : True { };
template<class F> struct IsFloat<Approximation<F>> : IsFloat<F> { };
template<class F> struct IsFloat<LowerBound<F>> : IsFloat<F> { };
template<class F> struct IsFloat<UpperBound<F>> : IsFloat<F> { };
template<class F> struct IsFloat<Bounds<F>> : IsFloat<F> { };
template<class F, class FE> struct IsFloat<Ball<F,FE>> : IsFloat<F> { };
template<class F> struct IsFloat<Value<F>> : IsFloat<F> { };
template<class F> struct IsFloat<Error<F>> : IsFloat<F> { };
template<> struct IsFloat<Dbl> : False { };

template<class F> concept AFloat = IsFloat<F>::value;

template<class T> struct IsNumber;
template<> struct IsNumber<Dbl> : True { };
template<> struct IsNumber<FloatDP> : True { };
template<> struct IsNumber<FloatMP> : True { };
template<class F> struct IsNumber<Approximation<F>> : IsNumber<F> { };
template<class F> struct IsNumber<Bounds<F>> : IsNumber<F> { };
template<class F, class FE> struct IsNumber<Ball<F,FE>> : IsNumber<F> { };
template<class F> struct IsNumber<Value<F>> : IsNumber<F> { };

template<class T> struct NumericTraits;

template<class T> struct NumericTraits<Positive<T>> : public NumericTraits<T> {
    typedef Positive<typename NumericTraits<T>::GenericType> GenericType;
    typedef Positive<typename NumericTraits<T>::OppositeType> OppositeType;
};

} // namespace Ariadne

#endif
