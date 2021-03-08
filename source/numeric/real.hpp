/***************************************************************************
 *            numeric/real.hpp
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

/*! \file numeric/real.hpp
 *  \brief
 */



#ifndef ARIADNE_REAL_HPP
#define ARIADNE_REAL_HPP

#include <functional>

#include "utility/typedefs.hpp"
#include "utility/pointer.hpp"
#include "utility/handle.hpp"

#include "numeric/logical.decl.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/float.decl.hpp"

#include "paradigm.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/sequence.hpp"

#include "dyadic.hpp"

namespace Ariadne {

//! A group
template<class G> class Group {
  public:
    G ineg() const; //!< \brief Inplace negation.
    friend G add(G,G); //!< \brief Addition
    friend G neg(G); //!< \brief Negation \relates Group
    friend OutputStream& operator<<(OutputStream& os, Group<G> const&); //!< \brief Write.
};

//! Symmetries
class Symmetry : public Group<Symmetry> {
  public:
    friend Symmetry sub(Symmetry,Symmetry); //!< \brief Subtraction
};

class Accuracy;

class Real;
class LowerReal;
class UpperReal;
class NaiveReal;
class PositiveReal;
class PositiveLowerReal;
class PositiveUpperReal;
class PositiveNaiveReal;
template<> struct IsNumber<Real> : True { };

class ValidatedReal;
class ValidatedUpperReal;
class ValidatedLowerReal;
class ApproximateReal;

class Dyadic;
template<class F> class Bounds;
template<class F> class LowerBound;
template<class F> class UpperBound;
template<class F> class Approximation;

using DyadicBounds = Bounds<Dyadic>;
using DyadicLowerBound = LowerBound<Dyadic>;
using DyadicUpperBound = UpperBound<Dyadic>;
using DyadicApproximation = Approximation<Dyadic>;

using PositiveDyadicBounds = PositiveBounds<Dyadic>;
using PositiveDyadicLowerBound = PositiveLowerBound<Dyadic>;
using PositiveDyadicUpperBound = PositiveUpperBound<Dyadic>;
using PositiveDyadicApproximation = PositiveApproximation<Dyadic>;


//! \ingroup NumericModule

extern const Real pi;
extern const Real infinity;

class RealInterface;
using LowerRealInterface = RealInterface;
using UpperRealInterface = RealInterface;
using NaiveRealInterface = RealInterface;

//! \ingroup NumericModule
//! \brief %Real number type \f$\R\f$ supporting elementary functions, comparisons and limits.
//! \details Effectively computable real numbers can be computed to arbitrary given accuracy as a (dyadic) rational approximation;
//! equivalently, as arbitrarily tight (dyadic) rational intervals.
//! They can be given as fast-converging Cauchy sequence of (dyadic) rationals,
//! or as a nested sequence of nested (dyadic) rational intervals whose intersection is the number itself.
//!
//! \b Example
//! The follow code creates a real number, computes it to an accuracy of 128 binary digits, a
//! and outputs the result as floating-point bounds with 96 bits of precision:
//! \snippet tutorials/numeric_usage.cpp Real_usage
//!
//! \sa LowerReal, UpperReal, NaiveReal
class Real
    : public Handle<const RealInterface>
//    , public DeclareRealOperations<Real,PositiveReal>
//    , public DeclareAnalyticFieldOperations<Real>
//    , public DeclareLatticeOperations<Real,PositiveReal>
//    , public DeclareComparisonOperations<Real,Kleenean,NegatedSierpinskian>
//    , public DefineFieldOperators<Real>
{
  public:
    typedef RealInterface Interface;
    typedef EffectiveTag Paradigm;
    typedef Real NumericType;
  public:

    //! \name Constructors
    //!@{
    Real(); //!< Default constructor yields the integer \c 0 as a real number.
    explicit Real(SharedPointer<const Interface>); //!< Construct from any class implementing the real number interface.
    explicit Real(ConvergentSequence<DyadicBounds> const&); //!< Construct from a sequence of dyadic bounds converging to a singleton intersection.
    explicit Real(FastCauchySequence<Dyadic> const&); //!< Construct from a fast convergent sequence of dyadic numbers i.e. \f$|w_m-w_n|\leq 2 ^\min(m,n)\f$.
    //!@}

    //! \name Conversion operators
    //!@{
    explicit Real(double) = delete; //!< Unsafe construction from a double is DEPRECATED

#ifdef DOXYGEN
    Real(Int n); //!< Convert from a builtin integer.
#else
    template<BuiltinUnsignedIntegral M> Real(M m);
    template<BuiltinSignedIntegral N> Real(N n);
#endif
    Real(ExactDouble d); //!< Construct from a double-precision value representing a number exactly.
    Real(Integer const& n); //!< Construct from an integer.
    Real(Dyadic const& d); //!< Construct from a dyadic number \a d=p/2<sup>q</sup> for integers \a p, \a q.
    Real(Decimal const& d); //!< Construct from a decimal number, given by its decimal expansion.
    Real(Rational const& q); //!< Construct from a rational number.

    explicit Real(FloatDPValue x); //!< DEPRECATED
    explicit Real(EffectiveNumber r); //!< DEPRECATED

    operator EffectiveNumber() const; //!< Convert to an effective number for computations.
    //!@}

    //! \name Convert to real number types describing less information.
    //!@{
    UpperReal upper() const; //!< A real number allowing only computation of upper bounds.
    LowerReal lower() const; //!< A real number allowing only computation of lower bounds.
    //FloatDPApproximation approx() const; //!< A real number allowing only computation of approximations with no guarantees on the accuracy.
    //!@}

    //! \name Computation of rigorous validated bounds on the number.
    //!@{
    //! Compute a concrete approximation with an error of at most \a 2<sup>-acc</sup> i.e. \a acc binary digits.
    ValidatedReal compute(Accuracy acc) const;
    //! Compute a concrete approximation with a bound on the error. The error bound converges to \a 0 as \a eff approaches infinity.
    //! The time taken should be roughly polynomial in \a eff.
    ValidatedReal compute(Effort eff) const;
    //! Compute a concrete approximation using dyadic bounds.
    DyadicBounds compute_get(Effort eff) const;
    //! Compute a concrete approximation using double-precision.
    FloatDPBounds compute_get(Effort eff, DoublePrecision pr) const;
    //! Compute a concrete approximation using the given precision.
    FloatMPBounds compute_get(Effort eff, MultiplePrecision pr) const;
    //! Compute a concrete approximation using double-precision.
    FloatDPBounds get(DoublePrecision pr) const;
    //! Compute a concrete approximation using the given precision.
    FloatMPBounds get(MultiplePrecision pr) const;
    //!@}

    //! \name Standard arithmetic operators
    //!@{
    friend Real operator+(Real const& r); //!< Unary plus.
    friend Real operator-(Real const& r); //!< Unary minus.
    friend Real operator+(Real const& r1, Real const& r2); //!< Plus.
    friend Real operator-(Real const& r1, Real const& r2); //!< Minus.
    friend Real operator*(Real const& r1, Real const& r2); //!< Times.
    friend Real operator/(Real const& r1, Real const& r2); //!< Divides.
    friend Real& operator+=(Real& r1, Real const& r2); //!< Inplace plus.
    friend Real& operator-=(Real& r1, Real const& r2); //!< Inplace minus.
    friend Real& operator*=(Real& r1, Real const& r2); //!< Inplace times.
    friend Real& operator/=(Real& r1, Real const& r2); //!< Inplace divides.
    //!@}

    //! \name Named arithmetical functions
    //!@{
    friend Real nul(Real const& r); //!< Constant \a 0.
    friend Real pos(Real const& r); //!< Identity \a +r.
    friend Real neg(Real const& r); //!< Negative \a -r.
    friend Real hlf(Real const& r); //!< Half \a r÷2.
    friend Real sqr(Real const& r); //!< Square \a r<sup>2</sup>.
    friend Real rec(Real const& r); //!< Reciprocal \a 1/r.
    friend Real add(Real const& r1, Real const& r2); //!< \brief Sum \a r1+r2.
    friend Real sub(Real const& r1, Real const& r2); //!< \brief Difference \a r1-r2.
    friend Real mul(Real const& r1, Real const& r2); //!< \brief Product \a r1×r2.
    friend Real div(Real const& r1, Real const& r2); //!< \brief Quotient \a r1÷r2.
    friend Real fma(Real const& r1, Real const& r2, Real const& r3); //!< \brief Fused multiply-and-add \a r1×r2+r3.
    friend Real pow(Real const& r, Int n); //!< \brief Power \a r<sup>n</sup>.
    //!@}

    //! \name Algebraic and transcendental functions
    //!@{
    friend Real sqrt(Real const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend Real exp(Real const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend Real log(Real const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend Real sin(Real const& r); //!< The sine of \a r.
    friend Real cos(Real const& r); //!< The cosine of \a r.
    friend Real tan(Real const& r); //!< The tangent of \a r, sin(\a r)/cos(\a r).
    friend Real asin(Real const& r); //!< The arc-sine of \a r.
    friend Real acos(Real const& r); //!< The arc-cosine of \a r.
    friend Real atan(Real const& r); //!< The arc-tangent of \a r.
    //!@}

    //! \name Lattice operations
    //!@{
    friend PositiveReal abs(Real const&); //!< Absolute value \a |r|.
    friend Real min(Real const& r1, Real const& r2); //!< The mimimum of \a r1 and \a r2.
    friend Real max(Real const& r1, Real const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //! \name Operations based on the metric structure.
    //!@{
    friend PositiveReal dist(Real const& r1, Real const& r2);
        //!< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    friend PositiveUpperReal mag(Real const& r);
        //!< An over-approximation to the absolute value of \a r.
    friend FloatDPError mag(Real const&, DoublePrecision);
    //!@}

    //! \name Comparison operations and operators.
    //!@{
    friend Kleenean lt(Real const& r1, Real const& r2); //!< Returns \c true if \a r1<r2, \c false if \a r1\>r2, and \c indetermiate if \a r1==r2.
    friend NegatedSierpinskian operator==(Real const& r1, Real const& r2); //!< Equality is undecidable and may only robustly be falsified.
    friend Sierpinskian operator!=(Real const& r1, Real const& r2); //!< Inequality is undecidable and may only robustly be verified.
//    friend Kleenean operator<=(Real const& r1, Real const& r2); //!< Comparison \c leq.
//    friend Kleenean operator>=(Real const& r1, Real const& r2); //!< Comparison .
    friend Kleenean operator< (Real const& r1, Real const& r2); //!< Comparison \c lt.
    friend Kleenean operator> (Real const& r1, Real const& r2); //!< Comparison .

    friend Kleenean sgn(Real const& r); //!< Returns \c true if \a r>0, \c false if \a r<0, and \c indeterminate if \a r==0.
    ValidatedKleenean check_sgn(Real r, Effort eff);
    friend Boolean operator>(Real const& r1, Pair<Real,Real> lu); //!< Given \a l<u, returns \a true if r>l and \a false if r<u.
    friend Boolean nondeterministic_greater(Real const& r, Rational const& a, Rational const& b); //!< Given \a a<b, returns \a true if r>a and \a false if r<b.

    friend Real choose(Case<LowerKleenean,Real> const& c1, Case<LowerKleenean,Real> const& c2);
        //!< A nonextensional choice, between the value of any of the valid cases.
    friend Real when(Case<UpperKleenean,Real> const& c1, Case<UpperKleenean,Real> const& c2);
        //!< A value equal to that of each of the valid cases.
    //!@}

    //! \name Limit operations.
    //!@{
    friend Real limit(ConvergentSequence<DyadicBounds> const& qbs);
        //!< Create a real number from a sequence of dyadic bounds whose width converges to 0.
    friend Real limit(FastCauchySequence<Dyadic> const& qs);
        //!< Create a real number from a sequence of dyadic numbers \f$q_n\f$ for which
        //!< \f$|q_{n_1}-q_{n_2}|\leq 2^{\min(n_1,n_2)}\f$
    friend Real limit(FastCauchySequence<Real> const& rs);
        //!< The limit of a sequence of real numbers \f$r_n\f$ for which
        //!< \f$|r_{n_1}-r_{n_2}|\leq 2^{\min(n_1,n_2)}\f$
    //!@}

    //! \name Rounding operations.
    //!@{
    friend Integer round(Real const& r);
        //!< Round to a "nearby" integer. This integer <em>need not</em> be the <em>closest</em> integer.
    //!@}


    //! \name Operations on the representation.
    //!@{
    friend Bool same(Real const&, Real const&); //!< Test equivalence of representation.
    //!@}

    //! \name Input/output operations
    //!@{
    friend OutputStream& operator<<(OutputStream& os, Real const& r); //!< Write to an output stream.
    //!@}
    double get_d() const;
  private:
    Real(std::int64_t n, Void*);
    Real(std::uint64_t m, Void*);
};

template<BuiltinUnsignedIntegral M> inline Real::Real(M m) : Real(std::uint64_t(m),nullptr) { }
template<BuiltinSignedIntegral N> inline Real::Real(N n) : Real(std::int64_t(n),nullptr) { }

Real choose(Case<LowerKleenean,Real> const& c1, Case<LowerKleenean,Real> const& c2);
Real when(Case<UpperKleenean,Real> const& c1, Case<UpperKleenean,Real> const& c2);

//! \ingroup NumericModule
//! \brief Lower real number type \f$\R_<\f$.
//! \details An effectively computable <em>lower real</em> is a real number for which it is possible to compute arbitrarily accurate (dyadic) rational lower bounds;
//! equivalently, a bounded increasing sequence of (dyadic) rationals converging to the number.
//! However, the convergence rate is not known, and there may be arbitrarily large upward jumps in the sequence.
//!
//! Multiplication and division are unsupported, since given \f$x \geq \underline{x}\f$ and \f$y\geq\underline{y}\f$,
//! we cannot deduce \f$x \times y\geq\underline{x}\times\underline{y}\f$ if the bounds are negative.
//! However, multiplication and reciprocation of positive lower reals <em>is</em> supported.
//!
//! %Positive lower real numbers are a natural type for representing the measure of an open set,
//! or the separation between two points in a metric space.
//! \sa Real, UpperReal, NaiveReal
class LowerReal
    : public Handle<const LowerRealInterface>
{
  public:
    typedef EffectiveTag Paradigm;
    typedef LowerReal NumericType;
  public:
    LowerReal(Real);
    explicit LowerReal(SharedPointer<const Interface>);
  public:
    //!@{
    //! \name Computation of rigorous validated lower bounds on the number
    ValidatedLowerReal compute(Effort eff) const;
    DyadicLowerBound compute_get(Effort eff) const;
    FloatDPLowerBound compute_get(Effort eff, DoublePrecision pr) const;
    FloatMPLowerBound compute_get(Effort eff, MultiplePrecision pr) const;
    //!@}
  public:
    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(LowerReal const& r) = delete; //!< \em No absolute value operator!
    friend LowerReal min(LowerReal const& r1, LowerReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend LowerReal max(LowerReal const& r1, LowerReal const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend LowerReal pos(LowerReal const& r); //!< Identity \a +r.
    friend UpperReal neg(LowerReal const& r); //!< Negative \a -r.
    friend LowerReal neg(UpperReal const& r); //!< Negative \a -r.
    friend LowerReal hlf(LowerReal const& r); //!< Half \a r÷2.
    friend LowerReal add(LowerReal const& r1, LowerReal const& r2); //!< \brief Sum \a r1+r2.
    friend LowerReal sub(LowerReal const& r1, UpperReal const& r2); //!< \brief Difference \a r1-r2.
    friend UpperReal sub(UpperReal const& r1, LowerReal const& r2); //!< \brief Difference \a r1-r2.

    friend NaiveReal mul(LowerReal const& r1, LowerReal const& r2) = delete; //!< \brief \em No multiplication operator, since non-monotone!
    friend NaiveReal div(LowerReal const& r1, UpperReal const& r2) = delete; //!< \brief \em No division operator, since non-monotone!
    friend NaiveReal rec(LowerReal const& r) = delete; //!< \brief \em No reciprocal operation, since non-monotone!
    //!@}

    //!@{
    //! \name Named arithmetical functions on positive numbers
    friend LowerReal mul(LowerReal const& r1, PositiveReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend LowerReal mul(PositiveReal const& r1, LowerReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend LowerReal div(LowerReal const& r1, PositiveReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend LowerReal div(PositiveReal const& r1, UpperReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend UpperReal div(PositiveReal const& r1, LowerReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.

    friend PositiveLowerReal add(PositiveLowerReal const& r1, PositiveLowerReal const& r2); //!< Addition \a r1+r2 preserves positivity.
    friend PositiveLowerReal mul(PositiveLowerReal const& r1, PositiveLowerReal const& r2); //!< Multiplication \a r1×r2 preserves positivity.
    friend PositiveLowerReal div(PositiveLowerReal const& r1, PositiveUpperReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveUpperReal div(PositiveUpperReal const& r1, PositiveLowerReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveLowerReal rec(PositiveUpperReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    friend PositiveUpperReal rec(PositiveLowerReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend PositiveLowerReal sqrt(PositiveLowerReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend PositiveLowerReal exp(LowerReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend LowerReal log(PositiveLowerReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend LowerReal atan(LowerReal const& r); //!< The arc-tangent of \a r.
    friend PositiveLowerReal atan(PositiveLowerReal const& r); //!< Arc-tangent preserves positivity

    friend NaiveReal sin(LowerReal const& r); //!< No sine operation, since non-monotone!
    friend NaiveReal cos(LowerReal const& r); //!< No cosine operation, since non-monotone!
    friend NaiveReal tan(LowerReal const& r); //!< No tangent operation, since non-monotone!
    //!@}

    //!@{
    //! \name Limit operations.
    friend LowerReal limit(IncreasingSequence<Dyadic> const& qlbs);
        //!< Create a real number from an increasing sequence of dyadic lower bounds.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, LowerReal const& r); //< Write to an output stream.
    //!@}

    //!@{
    //! \name Mixed operations
    friend Real min(LowerReal const& lr1, Real const& r2);
    friend Real min(Real const& r1, LowerReal const& lr2);
    //!@}
};

//! \ingroup NumericModule
//! \brief Upper real number type \f$\R_>\f$.
//! \details An effectively computable <em>upper real</em> is a real number for which it is possible to compute arbitrarily accurate (dyadic) rational upper bounds;
//! equivalently, a bounded decreasing sequence of (dyadic) rationals converging to the number.
//! However, the convergence rate is not known, and there may be arbitrarily large downward jumps in the sequence.
//!
//! Multiplication and division are unsupported, since given \f$x \leq \overline{x}\f$ and \f$y\leq\overline{y}\f$,
//! we cannot deduce \f$x \times y\leq\overline{x}\times\overline{y}\f$.
//! However, multiplication and reciprocation of positive upper reals <em>is</em> supported.
//!
//! %Positive upper real numbers are a natural type for representing upper bounds of a distance or norm,
//! including the radius of a ball, or the measure of an closed set.
//! \sa Real, LowerReal, NaiveReal
class UpperReal
    : public Handle<const UpperRealInterface>
{
  public:
    typedef EffectiveTag Paradigm;
    typedef UpperReal NumericType;
  public:
    UpperReal(Real);
    explicit UpperReal(SharedPointer<const Interface>);
  public:
    //!@{
    //! \name Computation of rigorous validated upper bounds on the number
    ValidatedUpperReal compute(Effort eff) const;
    DyadicUpperBound compute_get(Effort eff) const;
    FloatDPUpperBound compute_get(Effort eff, DoublePrecision pr) const;
    FloatMPUpperBound compute_get(Effort eff, MultiplePrecision pr) const;
    //!@}

  public:
    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(UpperReal const& r) = delete; //!< \em No absolute value operator!
    friend UpperReal min(UpperReal const& r1, UpperReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend UpperReal max(UpperReal const& r1, UpperReal const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend UpperReal pos(UpperReal const& r); //!< Identity \a +r.
    friend LowerReal neg(UpperReal const& r); //!< Negative \a -r.
    friend UpperReal neg(LowerReal const& r); //!< Negative \a -r.
    friend UpperReal hlf(UpperReal const& r); //!< Half \a r÷2.
    friend UpperReal add(UpperReal const& r1, UpperReal const& r2); //!< \brief Sum \a r1+r2.
    friend UpperReal sub(UpperReal const& r1, LowerReal const& r2); //!< \brief Difference \a r1-r2.
    friend LowerReal sub(LowerReal const& r1, UpperReal const& r2); //!< \brief Difference \a r1-r2.

    friend NaiveReal mul(UpperReal const& r1, UpperReal const& r2) = delete; //!< \brief \em No multiplication operator, since non-monotone!
    friend NaiveReal mul(UpperReal const& r1, LowerReal const& r2) = delete; //!< \brief \em No division operator, since non-monotone!
    friend NaiveReal rec(UpperReal const& r) = delete; //!< \brief \em No reciprocal operation, since non-monotone!
    //!@}

    //!@{
    //! \name Named arithmetical functions on positive numbers
    friend UpperReal mul(UpperReal const& r1, PositiveReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend UpperReal mul(PositiveReal const& r1, UpperReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend UpperReal div(UpperReal const& r1, PositiveReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend UpperReal div(PositiveReal const& r1, LowerReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend LowerReal div(PositiveReal const& r1, UpperReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.

    friend PositiveUpperReal add(PositiveUpperReal const& r1, PositiveUpperReal const& r2); //!< Addition \a r1+r2 preserves positivity.
    friend PositiveUpperReal mul(PositiveUpperReal const& r1, PositiveUpperReal const& r2); //!< Multiplication \a r1×r2 preserves positivity.
    friend PositiveUpperReal div(PositiveUpperReal const& r1, PositiveLowerReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveLowerReal div(PositiveLowerReal const& r1, PositiveUpperReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveUpperReal rec(PositiveLowerReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    friend PositiveLowerReal rec(PositiveUpperReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend PositiveUpperReal sqrt(PositiveUpperReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend PositiveUpperReal exp(UpperReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend UpperReal log(PositiveUpperReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend UpperReal atan(UpperReal const& r); //!< The arc-tangent of \a r.
    friend PositiveUpperReal atan(PositiveUpperReal const& r); //!< Arc-tangent preserves positivity
    friend NaiveReal sin(UpperReal const& r); //!< No sine operation, since non-monotone!
    friend NaiveReal cos(UpperReal const& r); //!< No cosine operation, since non-monotone!
    friend NaiveReal tan(UpperReal const& r); //!< No tangent operation, since non-monotone!
    //!@}

    //!@{
    //! \name Limit operations.
    friend UpperReal limit(DecreasingSequence<Dyadic> const& qubs);
        //!< Create a real number from an decreasing sequence of dyadic upper bounds.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, UpperReal const& r); //< Write to an output stream.
    //!@}

    //!@{
    //! \name Mixed operations
    friend Real max(UpperReal const& ur1, Real const& r2);
    friend Real max(Real r1, UpperReal const& ur2);
    //!@}
};

//! \ingroup NumericModule
//! \brief %Real number type defined as limits of convergent sequences of (dyadic) rationals, without bounds on the convergence rate.
//! \details It is possible to compute arbitrarily accurate (dyadic) rational approximations, but the convergence rate is not known,
//! so at no point in the the computation can anything concrete be deduced about the value of the number.
//!
//! In principle useless for rigorous computation, but quickly-computed approximations to real numbers may be useful for preconditioning rigorous algorithms.
//! \sa Real, LowerReal, UpperReal
class NaiveReal
    : public Handle<const NaiveRealInterface>
    , public DeclareRealOperations<NaiveReal,PositiveNaiveReal>
    , public DeclareAnalyticFieldOperations<NaiveReal>
    , public DeclareLatticeOperations<NaiveReal,PositiveNaiveReal>
    , public DeclareComparisonOperations<NaiveReal,ApproximateKleenean>
    , public DefineFieldOperators<NaiveReal>
{
  public:
    typedef NaiveRealInterface Interface;
    typedef ApproximateTag Paradigm;
    typedef NaiveReal NumericType;
  public:
    //!@{
    //! \name Constructors
    NaiveReal(); //!< Default constructor yields the integer \c 0 as a real number.
    explicit NaiveReal(SharedPointer<const Interface>); //!< Construct from any class implementing the real number interface.
    //explicit NaiveReal(ConvergentSequence<DyadicBounds> const&); //!< Construct from a sequence of dyadic bounds converging to a singleton intersection.
    //!@}

    //!@{
    //! \name Conversion operators
    NaiveReal(Dbl d); //!< Construction from a builtin double-precision approximation.
    NaiveReal(Integer const& n); //!< Construct from an integer.
    NaiveReal(Dyadic const& d); //!< Construct from a dyadic number \a d=p/2<sup>q</sup> for integers \a p, \a q.
    NaiveReal(Decimal const& d); //!< Construct from a decimal number, given by its decimal expansion.
    NaiveReal(Rational const& q); //!< Construct from a rational number.

    NaiveReal(Real const& r); //!< Construct from an (effective) real number.
    NaiveReal(LowerReal const& lr); //!< Construct from a lower real number.
    NaiveReal(UpperReal const& ur); //!< Construct from an upper real number.

//    operator Number<Tag>() const; //!< Convert to an effective number for computations.
    //!@}

    //!@{
    //! \name Computation of approximations to the number.
    //! Compute an approximation with no guarantees on the error. The error converges to \a 0 as \a eff approaches infinity.
    //! The time taken should be roughly polynomial in \a eff.
    ApproximateReal compute(Effort eff) const;
    //!@}

    //!@{
    //! \name Standard arithmetic operators
    friend NaiveReal operator+(NaiveReal const& r); //!< Unary plus.
    friend NaiveReal operator-(NaiveReal const& r); //!< Unary minus.
    friend NaiveReal operator+(NaiveReal const& r1, NaiveReal const& r2); //!< Plus.
    friend NaiveReal operator-(NaiveReal const& r1, NaiveReal const& r2); //!< Minus.
    friend NaiveReal operator*(NaiveReal const& r1, NaiveReal const& r2); //!< Times.
    friend NaiveReal operator/(NaiveReal const& r1, NaiveReal const& r2); //!< Divides.
    friend NaiveReal& operator+=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace plus.
    friend NaiveReal& operator-=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace minus.
    friend NaiveReal& operator*=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace times.
    friend NaiveReal& operator/=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace divides.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend NaiveReal pos(NaiveReal const& r); //!< Identity \a +r.
    friend NaiveReal neg(NaiveReal const& r); //!< Negative \a -r.
    friend NaiveReal hlf(NaiveReal const& r); //!< Half \a r÷2.
    friend NaiveReal sqr(NaiveReal const& r); //!< Square \a r<sup>2</sup>.
    friend NaiveReal rec(NaiveReal const& r); //!< Reciprocal \a 1/r.
    friend NaiveReal add(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Sum \a r1+r2.
    friend NaiveReal sub(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Difference \a r1-r2.
    friend NaiveReal mul(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Product \a r1×r2.
    friend NaiveReal div(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Quotient \a r1÷r2.
    friend NaiveReal fma(NaiveReal const& r1, NaiveReal const& r2, NaiveReal const& r3); //!< \brief Fused multiply-and-add \a r1×r2+r3.
    friend NaiveReal pow(NaiveReal const& r, Int n); //!< \brief Power \a r<sup>n</sup>.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend NaiveReal sqrt(NaiveReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend NaiveReal exp(NaiveReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend NaiveReal log(NaiveReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend NaiveReal sin(NaiveReal const& r); //!< The sine of \a r.
    friend NaiveReal cos(NaiveReal const& r); //!< The cosine of \a r.
    friend NaiveReal tan(NaiveReal const& r); //!< The tangent of \a r, sin(\a r)/cos(\a r).
    friend NaiveReal atan(NaiveReal const& r); //!< The arc-tangent of \a r.
    //!@}

    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(NaiveReal const&); //!< Absolute value \a |r|.
    friend NaiveReal min(NaiveReal const& r1, NaiveReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend NaiveReal max(NaiveReal const& r1, NaiveReal const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! Operations based on the metric structure.
    friend PositiveNaiveReal dist(NaiveReal const& r1, NaiveReal const& r2);
        //< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    //!@}

    //!@{
    //! \name Comparison operations and operators.
    friend ApproximateKleenean leq(NaiveReal const& r1, NaiveReal const& r2);
        //!< Returns \c likely if \a r1<r2, \c unlikely if \a r1\>r2, and \c indetermiate if \a r1==r2.
    //!@}

    //!@{
    //! \name Limit operations.
    friend NaiveReal limit(ConvergentSequence<Dyadic> const& qas);
        //!< Create a real number from an converging sequence of dyadic approximants.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, NaiveReal const& r); //< Write to an output stream.
    //!@}
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveReal : public Real
{
  public:
    PositiveReal() : Real() { }
    template<ConvertibleTo<Real> Y>
        PositiveReal(Positive<Y> const& p) : Real(p) { }
    explicit PositiveReal(Real const& r) : Real(r) { }
    PositiveBounds<Dyadic> compute_get(Effort) const;
  public:
    PositiveReal max(PositiveReal const&, PositiveReal const&);
    PositiveReal max(Real const&, PositiveReal const&);
    PositiveReal max(PositiveReal const&, Real const&);

    PositiveReal min(PositiveReal const&, PositiveReal const&);

    PositiveReal rec(PositiveReal const&);
    PositiveReal add(PositiveReal const&, PositiveReal const&);
    PositiveReal mul(PositiveReal const&, PositiveReal const&);
    PositiveReal div(PositiveReal const&, PositiveReal const&);
  public:
    friend PositiveReal min(PositiveReal const& pr1, PositiveReal const& pr2);
    friend PositiveReal add(PositiveReal const& pr1, PositiveReal const& pr2);
    friend PositiveReal mul(PositiveReal const& pr1, PositiveReal const& pr2);
    friend PositiveReal div(PositiveReal const& pr1, PositiveReal const& pr2);
    friend PositiveReal rec(PositiveReal const& pr);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveLowerReal : public LowerReal
{
  public:
    PositiveLowerReal(PositiveReal r) : LowerReal(r) { }
    explicit PositiveLowerReal(LowerReal r) : LowerReal(r) { }
    PositiveLowerBound<Dyadic> compute_get(Effort) const;
  public:
    PositiveLowerReal rec(PositiveUpperReal const&);
    PositiveUpperReal rec(PositiveLowerReal const&);
    PositiveLowerReal add(PositiveLowerReal const&, PositiveLowerReal const&);
    PositiveLowerReal mul(PositiveLowerReal const&, PositiveLowerReal const&);
    PositiveLowerReal div(PositiveLowerReal const&, PositiveUpperReal const&);
    PositiveUpperReal div(PositiveUpperReal const&, PositiveLowerReal const&);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveUpperReal : public UpperReal
{
  public:
    PositiveUpperReal(PositiveReal r) : UpperReal(r) { }
    explicit PositiveUpperReal(UpperReal r) : UpperReal(r) { }
    PositiveUpperBound<Dyadic> compute_get(Effort) const;
  public:
    PositiveUpperReal rec(PositiveLowerReal const&);
    PositiveLowerReal rec(PositiveUpperReal const&);
    PositiveUpperReal add(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal mul(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal div(PositiveUpperReal const&, PositiveLowerReal const&);
    PositiveLowerReal div(PositiveLowerReal const&, PositiveUpperReal const&);
};

PositiveUpperReal rec(PositiveLowerReal plr);
PositiveLowerReal rec(PositiveUpperReal pur);
PositiveLowerReal add(PositiveLowerReal plr1, PositiveLowerReal plr2);
PositiveUpperReal add(PositiveUpperReal pur1, PositiveUpperReal pur2);
PositiveLowerReal mul(PositiveLowerReal plr1, PositiveLowerReal plr2);
PositiveUpperReal mul(PositiveUpperReal pur1, PositiveUpperReal pur2);
PositiveLowerReal div(PositiveLowerReal plr1, PositiveUpperReal pur2);
PositiveUpperReal div(PositiveUpperReal pur1, PositiveLowerReal plr2);

PositiveReal cast_positive(Real const& x);

//LowerReal add(LowerReal const& lr1, UpperReal const& ur2);
//UpperReal add(UpperReal const& ur1, LowerReal const& lr2);

ValidatedKleenean check_sgn(Real r, Effort eff);

class ValidatedRealInterface;

//! \ingroup NumericModule
//! \brief A generic class representing rigorous bounds on a real number.
//! \see Real
class ValidatedReal
    : public Handle<ValidatedRealInterface>
{
  public:
    typedef ValidatedRealInterface Interface;

    ValidatedReal(DyadicBounds const&);
    operator DyadicBounds() const;

    Dyadic value();
    Dyadic error();
    Dyadic lower();
    Dyadic upper();

    //! \brief Get dyadic bounds for the number.
    //! \details It is not always clear how this function can be implemented,
    //!   and it may be removed or modified in future versions.
    DyadicBounds get() const;
    //! \brief Get the bounds for the number, representing in double precision.
    FloatDPBounds get(DoublePrecision) const;
    //! \brief Get the bounds for the number, representing to precision \a pr.
    //! Note that increasing the accuracy typically does not yield arbitrarily
    //! tight bounds, as the object is already an approximation to finite accuracy.
    FloatMPBounds get(MultiplePrecision pr) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream&, ValidatedReal const&);

    friend ValidatedKleenean sgn(DyadicBounds const& r);
    friend ValidatedKleenean operator==(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator!=(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator<=(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator>=(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator< (DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator> (DyadicBounds const& vr1, DyadicBounds const& vr2);
};

class ValidatedLowerReal
    : public Handle<ValidatedRealInterface>
{
  public:
    ValidatedLowerReal(DyadicLowerBound const&);
    ValidatedLowerReal(ValidatedReal const&);
    DyadicLowerBound get() const;
    FloatDPLowerBound get(DoublePrecision) const;
    FloatMPLowerBound get(MultiplePrecision) const;
    friend OutputStream& operator<<(OutputStream&, ValidatedLowerReal const&);
};

class ValidatedUpperReal
    : public Handle<ValidatedRealInterface>
{
  public:
    ValidatedUpperReal(DyadicUpperBound const&);
    ValidatedUpperReal(ValidatedReal const&);
    DyadicUpperBound get() const;
    FloatDPUpperBound get(DoublePrecision) const;
    FloatMPUpperBound get(MultiplePrecision) const;
    friend OutputStream& operator<<(OutputStream&, ValidatedUpperReal const&);
};

class ApproximateRealInterface;

//! \ingroup NumericModule
//! \brief A generic class representing an approximation to a real number
//! of unknown accuracy.
//! \see Real, ValidatedReal
class ApproximateReal
    : public Handle<ApproximateRealInterface>
{
  public:
    typedef ApproximateRealInterface Interface;

    ApproximateReal(DyadicApproximation const&);
    //! \brief Get a dyadic approximation to the number.
    DyadicApproximation get() const;
    //! \brief Get the approximation to the number, converting to double-precision.
    FloatDPApproximation get(DoublePrecision) const;
    //! \brief Get the approximation to the number, converting to a representation
    //! with precision \a pr.
    //! Note that increasing the precision typically does not yield arbitrarily
    //! accurate approximation, as the object already has some error.
    FloatMPApproximation get(MultiplePrecision pr) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream&, ApproximateReal const&);
};


} // namespace Ariadne

#endif
