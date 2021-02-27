/***************************************************************************
 *            numeric/dyadic.hpp
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

/*! \file numeric/dyadic.hpp
 *  \brief Dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_HPP
#define ARIADNE_DYADIC_HPP

#include "numeric/gmp.hpp"

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "utility/handle.hpp"
#include "utility/writable.hpp"

#include "numeric/concepts.hpp"
#include "numeric/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/bounds.hpp"

namespace Ariadne {

class ExactDouble;
class Dyadic;

class InfinityException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

//! \ingroup NumericModule
//! \related FloatDP, ExactIntervalType
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
class Dyadic
//    : DeclareRingOperations<Dyadic>
//    , DeclareLatticeOperations<Dyadic,Dyadic>
//    , DeclareComparisonOperations<Dyadic,Boolean,Boolean>
//    , DefineRingOperators<Dyadic>
//    , DefineComparisonOperators<Dyadic,Boolean,Boolean>
    : ProvideOperators<Dyadic>
{
    static Writer<Dyadic> _default_writer;
  public:
//    template<class W> requires BaseOf<WriterInterface<Dyadic>,W> static Void set_default_writer(W w) {
//        _default_writer=std::make_shared<W>(std::move(w)); }
//    static Void set_default_writer(Writer<Dyadic> w) { _default_writer=w.managed_pointer(); }
    static Void set_default_writer(Writer<Dyadic> w) { _default_writer=w; }
    static Writer<Dyadic> default_writer() { return _default_writer; }
  public:
    mpf_t _mpf;
  public:
    typedef ExactTag Paradigm;

    //! \brief Construct a Dyadic number from a GNU mpf object.
    explicit Dyadic (mpf_t mpf);
    //! \brief Construct the Dyadic number \a p/2<sup>q</sup>.
    Dyadic (Integer const& p, Nat q);
    Dyadic (Integer const& p, Int q) = delete; // Can only construct a dyadic p/2^q for a positive (unsigned) value q, such at uint or Nat.
    Dyadic (Integer const& p, Natural q);
    //! \brief Destructor.
    ~Dyadic();
    //! \brief Default constructor creates the number 0 (zero).
    Dyadic();
    //! \brief Copy constructor.
    Dyadic(Dyadic const& n);
    Dyadic(Dyadic&& n);
    //! \brief Assignment constructor.
    Dyadic& operator=(Dyadic const& n);
    Dyadic& operator=(Dyadic&& n);
    //! \brief Convert from a built-in integer.
    template<BuiltinIntegral N> Dyadic(N n);
    //! \brief Convert from an exact double-precision number.
    Dyadic(const ExactDouble& d);
    //! \brief Convert from an integer.
    Dyadic(const Integer& z);
    //! \brief Convert from a power of two.
    Dyadic(const TwoExp& t);
    //! \brief Explicit construction from a built-in double-precision value.
    //! \details Tests to ensure that the number is not 'accidentally' created from a rounded version of a string literal,
    //! by comparing the input with it's single-precision approximation.
    explicit Dyadic(double x);
    //! \brief Convert to a generic number.
    operator ExactNumber () const;
    //! \brief A representation of ±∞ or NaN.
    static Dyadic inf(Sign sgn);
    //! \brief A representation of +∞.
    static Dyadic inf();
    //! \brief A representation of NaN (not-a-number).
    static Dyadic nan();

    //! \brief The smallest integer \a p such that \a x=p/2<sup>q</sup>
    Integer mantissa() const;
    //! \brief The (negative) integer \a -q such that \a x=p/2<sup>q</sup>
    Int exponent() const;
    //! \brief A double-precision approximateion.
    double get_d() const;
    mpf_t const& get_mpf() const;
    //! \brief Convert a floating-point literal to Dyadic i.e. long binary format.
    friend Dyadic operator"" _bin(long double x);
    //! \brief Convert a floating-point literal to Dyadic.
    friend Dyadic operator"" _dyadic(long double x);

    friend Dyadic nul(Dyadic const&);
    friend Dyadic pos(Dyadic const&);
    friend Dyadic neg(Dyadic const&);
    friend Dyadic sqr(Dyadic const& q);
    friend Dyadic add(Dyadic const&, Dyadic const&);
    friend Dyadic sub(Dyadic const&, Dyadic const&);
    friend Dyadic mul(Dyadic const&, Dyadic const&);
    //! \brief Halve the number.
    friend Dyadic hlf(Dyadic const&);
    //| \brief Power of a number (m always positive).
    friend Dyadic pow(Dyadic const&, Int);
    friend Dyadic pow(Dyadic const&, Nat);

    friend Rational rec(Rational const&);
    friend Rational div(Rational const&, Rational const&);
//    friend Rational operator/(Rational const&, Rational const&);

    friend Real sqrt(Real const&); //!< <p/>
    friend Real exp(Real const&); //!< <p/>
    friend Real log(Real const&); //!< <p/>
    friend Real sin(Real const&); //!< <p/>
    friend Real cos(Real const&); //!< <p/>
    friend Real tan(Real const&); //!< <p/>
    friend Real asin(Real const&); //!< <p/>
    friend Real acos(Real const&); //!< <p/>
    friend Real atan(Real const&); //!< <p/>

    friend Dyadic abs(Dyadic const&);
    friend Dyadic max(Dyadic const&, Dyadic const&);
    friend Dyadic min(Dyadic const&, Dyadic const&);

    friend Comparison cmp(Dyadic const&, Dyadic const&);

    //! \brief The sign of the number.
    friend Sign sgn(Dyadic const&);
    //! \brief Round down to the nearest lower integer.
    friend Integer floor(Dyadic const&);
    //! \brief Round to the nearest integer. Rounding of halves is implementation-dependent.
    friend Integer round(Dyadic const&);
    //! \brief Round up to the nearest higher integer.
    friend Integer ceil(Dyadic const&);

    //! \brief Tests whether the value is NaN (not-a-number).
    friend Bool is_nan(Dyadic const& w);
    //! \brief Tests whether the value is ±∞.
    friend Bool is_inf(Dyadic const& w);
    //! \brief Tests whether the value is finite.
    friend Bool is_finite(Dyadic const& w);
    //! \brief Tests whether the value is zero.
    friend Bool is_zero(Dyadic const& w);

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, Dyadic const& x);
};

class DecimalWriter : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os, Dyadic const& w) const final override;
};
class FractionWriter : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os, Dyadic const& w) const final override;
};
template<> class RepresentationWriter<Dyadic> : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os, Dyadic const& w) const final override;
};


template<BuiltinIntegral N> inline Dyadic::Dyadic(N n) : Dyadic(Integer(n)) { }


template<> class Positive<Dyadic> : public Dyadic {
  public:
    Positive<Dyadic>() : Dyadic() { }
    template<BuiltinUnsignedIntegral M> Positive<Dyadic>(M m) : Dyadic(m) { }
    Positive<Dyadic>(int n) = delete;
    explicit Positive<Dyadic>(Dyadic const& z) : Dyadic(z) { assert(z>=0); }
};
inline Positive<Dyadic> cast_positive(Dyadic const& w) { return Positive<Dyadic>(w); }

using PositiveDyadic = Positive<Dyadic>;

//! \relates Dyadic
//! \name Type synonyms
//!@{
using DyadicBall = Ball<Dyadic,Dyadic>; //!< Alias for ball about a number with dyadic value and error.
using DyadicBounds = Bounds<Dyadic>; //!< Alias for dyadic bounds on a number.
using DyadicUpperBound = UpperBound<Dyadic>; //!< Alias for dyadic upper bound for a number.
using DyadicLowerBound = LowerBound<Dyadic>; //!< Alias for dyadic lower bound for a number.
using DyadicApproximation = Approximation<Dyadic>; //!< Alias for dyadic approximateion to a number.

using PositiveDyadicApproximation = Positive<Approximation<Dyadic>>; //!< <p/>
using PositiveDyadicLowerBound = Positive<LowerBound<Dyadic>>; //!< <p/>
using PositiveDyadicUpperBound = Positive<UpperBound<Dyadic>>; //!< <p/>
using PositiveDyadicBounds = Positive<Bounds<Dyadic>>; //!< <p/>
//!@}

template<class F> class Bounds;
template<class F, class FE> class Ball;

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Bounds<Y>;
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Ball<Y>;

using DyadicBounds = Bounds<Dyadic>; //!< Alias for dyadic bounds on a number. //!< \ingroup NumericModule
using DyadicBall = Ball<Dyadic,Dyadic>;

static_assert(LatticeRing<Dyadic>);
static_assert(not RoundedRing<Dyadic>);

template class Bounds<Dyadic>;
template class Ball<Dyadic>;

using DyadicApproximation = Approximation<Dyadic>;

template<class F> class Bounds;
template<class F, class FE> class Ball;

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Bounds<Y>;
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Ball<Y>;

using DyadicBounds = Bounds<Dyadic>; //!< Alias for dyadic bounds on a number. //!< \ingroup NumericModule
using DyadicBall = Ball<Dyadic,Dyadic>;

static_assert(LatticeRing<Dyadic>);
static_assert(not RoundedRing<Dyadic>);

using DyadicApproximation = Approximation<Dyadic>;

inline Dyadic operator"" _dyadic(long double x) { return Dyadic(static_cast<double>(x)); }
inline Dyadic operator"" _dy(long double x) { return operator"" _dyadic(x); }
inline Dyadic operator"" _q2(long double x) { return operator"" _dyadic(x); }
inline Dyadic operator"" _bin(long double x) { return operator"" _dyadic(x); }

Comparison cmp(Dyadic const& x1, Dyadic const& x2);
Dyadic make_dyadic(unsigned long long int n);


} // namespace Ariadne

#endif
