/***************************************************************************
 *            numeric/rational.hpp
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

/*! \file numeric/rational.hpp
 *  \brief
 */



#ifndef ARIADNE_RATIONAL_HPP
#define ARIADNE_RATIONAL_HPP

#include "../numeric/gmp.hpp"
#include "../utility/typedefs.hpp"
#include "../utility/metaprogramming.hpp"
#include "../utility/string.hpp"
#include "../numeric/integer.hpp"
#include "../numeric/arithmetic.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

namespace Ariadne {

class Int64;
class FloatDP;

enum class Comparison : char;


/************ Rational *******************************************************/

//! \ingroup NumericModule
//! \brief %Rational numbers.
class Rational
    : DeclareFieldOperations<Rational>
    , DeclareLatticeOperations<Rational,Rational>
    , DeclareComparisonOperations<Rational,Boolean,Boolean>
    , DefineFieldOperators<Rational>
    , DefineComparisonOperators<Rational,Boolean,Boolean>
    , DeclareTranscendentalOperations<Real>
{
  public:
    mpq_t _mpq;
  public:
    typedef ExactTag Paradigm;
    typedef Rational NumericType;
  public:
    ~Rational();
    Rational();
    Rational(const Integer&, const Integer&);
    template<BuiltinIntegral N> Rational(N n);
    Rational(Int64);
    explicit Rational(FloatDP const&);
    Rational(const ExactDouble&);
    Rational(const Integer&);
    Rational(const Dyadic&);
    explicit Rational(const String&);
    explicit Rational(const FloatDPValue&);
    explicit Rational(const mpq_t);
    Rational(const Rational&);
    Rational(Rational&&);
    Rational& operator=(const Rational&);
    Rational& operator=(Rational&&);
    operator ExactNumber () const;

    static Rational inf(Sign sgn);
    static Rational inf();
    static Rational nan();

    Integer get_num() const;
    Integer get_den() const;
    Integer numerator() const;
    Natural denominator() const;
    friend Rational operator/(Integer const& z1, Integer const& z2);

    friend Real sqrt(Real const&);
    friend Real exp(Real const&);
    friend Real log(Real const&);
    friend Real sin(Real const&);
    friend Real cos(Real const&);
    friend Real tan(Real const&);
    friend Real asin(Real const&);
    friend Real acos(Real const&);
    friend Real atan(Real const&);

    friend Sign sgn(Rational const& q);
    friend Integer floor(Rational const&);
    friend Integer round(Rational const&);
    friend Integer ceil(Rational const&);

    friend Bool is_nan(Rational const& q);
    friend Bool is_inf(Rational const& q);
    friend Bool is_finite(Rational const& q);
    friend Bool is_zero(Rational const& q);

    friend Comparison cmp(Rational const& q1, Rational const& q2);
    friend Comparison cmp(Rational const& q1, ExactDouble const& d2);
    friend Comparison cmp(ExactDouble const& d1, Rational const& q2);

    friend OutputStream& operator<<(OutputStream& os, Rational const& q);
    friend InputStream& operator>>(InputStream& os, Rational& q);
    friend Rational operator"" _q(long double x);
  public:
    double get_d() const;
    mpq_t const& get_mpq() const;
  private:
    friend class Dyadic;
};
template<> struct IsNumericType<Rational> : True { };
Rational operator"" _q(unsigned long long int n);
Rational operator"" _q(long double x);

template<BuiltinIntegral N> inline Rational::Rational(N n) : Rational(Int64(n)) { }

OutputStream& write(OutputStream& os, mpz_t const z);
InputStream& operator>>(InputStream& is, Rational& q1);

using RationalBounds = Bounds<Rational>; //!< Alias for rational bounds on a number. //!< \ingroup NumericModule

template<> class Bounds<Rational> {
    Rational _l, _u;
  public:
    Bounds<Rational>(Rational q) : _l(q), _u(q) { }
    Bounds<Rational>(Rational l, Rational u) : _l(l), _u(u) { }
    template<class X> requires Constructible<Rational,X> Bounds<Rational>(Bounds<X> const& x)
        : Bounds<Rational>(Rational(x.lower_raw()),Rational(x.upper_raw())) { }
    operator ValidatedNumber() const;
    Bounds<Rational> pm(Rational e) { return RationalBounds(_l-e,_u+e); }
    Rational lower() const { return _l; }
    Rational upper() const { return _u; }
    Rational lower_raw() const { return _l; }
    Rational upper_raw() const { return _u; }
    friend Bounds<Rational> operator+(RationalBounds const& q) { return RationalBounds(+q._l,+q._u); }
    friend Bounds<Rational> operator-(RationalBounds const& q) { return RationalBounds(-q._u,-q._l); }
    friend Bounds<Rational> operator+(RationalBounds const& q1, RationalBounds const& q2) { return RationalBounds(q1._l+q2._l,q1._u+q2._u); }
    friend Bounds<Rational> operator-(RationalBounds const& q1, RationalBounds const& q2) { return RationalBounds(q1._l-q2._u,q1._u-q2._l); }
    friend Bounds<Rational> operator*(RationalBounds const& q1, RationalBounds const& q2);
    friend Bounds<Rational> operator/(RationalBounds const& q1, RationalBounds const& q2);
    friend Bounds<Rational> sqr(RationalBounds const& q);
    friend Bounds<Rational> rec(RationalBounds const& q);
    friend RationalBounds abs(RationalBounds const& q) { return RationalBounds(max(min(q._l,-q._u),0),max(-q._l,q._u)); }
    friend Bounds<Rational> max(RationalBounds const& q1, RationalBounds const& q2) { return RationalBounds(max(q1._l,q2._l),max(q1._u,q2._u)); }
    friend Bounds<Rational> min(RationalBounds const& q1, RationalBounds const& q2) { return RationalBounds(min(q1._l,q2._l),min(q1._u,q2._u)); }
    friend ValidatedKleenean operator<(RationalBounds const& q1, RationalBounds const& q2) {
        if (q1._u<q2._l) { return true; } else if (q1._l >= q2._u) { return false; } else { return indeterminate; } }
    friend Boolean refines(RationalBounds const& q1, RationalBounds const& q2) { return q1._l>=q2._l and q1._u<=q2._u; }
    friend RationalBounds refinement(RationalBounds const& q1, RationalBounds const& q2) { return RationalBounds(max(q1._l,q2._l),min(q1._u,q2._u)); }
    friend OutputStream& operator<<(OutputStream& os, RationalBounds y) { return os << "[" << y._l << ":" << y._u << "]"; }
};


} // namespace Ariadne

#endif
