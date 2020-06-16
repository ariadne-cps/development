/***************************************************************************
 *            numeric/integer.hpp
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

/*! \file numeric/integer.hpp
 *  \brief
 */



#ifndef ARIADNE_INTEGER_HPP
#define ARIADNE_INTEGER_HPP

#include "../numeric/gmp.hpp"

#include <cassert>

#include "../utility/typedefs.hpp"
#include "../utility/metaprogramming.hpp"
#include "../numeric/sign.hpp"
#include "../numeric/logical.hpp"
#include "../numeric/arithmetic.hpp"
#include "../numeric/number.decl.hpp"

namespace Ariadne {

/************  Ints ********************************************************/

class Nat32 {
    uint32_t _m;
  public:
    Nat32() : _m(0u) { }
    template<BuiltinUnsignedIntegral M> Nat32(M m) : _m(m) { assert(_m==m); }
    template<BuiltinSignedIntegral N> Nat32(N n) : _m(n) { assert(n>=0); assert((int64_t)_m==n); }
    uint32_t get_ui() const { return _m; }
};

class Nat64 {
    uint64_t _m;
  public:
    Nat64() : _m(0u) { }
    template<BuiltinUnsignedIntegral M> Nat64(M m) : _m(m) { assert(_m==m); }
    template<BuiltinSignedIntegral N> Nat64(N n) : _m(static_cast<uint64_t>(n)) { assert(n>=0); assert((int64_t)_m==n);
        assert(uint64_t(int64_t(_m))==_m); }
    uint64_t get_ui() const { return _m; }
};

class Int32 {
    int32_t _n;
  public:
    Int32() : _n(0) { }
    template<BuiltinUnsignedIntegral M> Int32(M m) : _n(static_cast<int32_t>(m)) { assert(_n>=0); assert((uint32_t)_n==m); }
    template<BuiltinSignedIntegral N> Int32(N n) : _n(n) { assert(_n==n); }
    int32_t get_si() const { return _n; }
};

class Int64 {
    int64_t _n;
  public:
    Int64() : _n(0) { }
    template<BuiltinUnsignedIntegral M> Int64(M m) : _n(static_cast<int64_t>(m)) { assert(_n>=0); assert((uint64_t)_n==m); }
    template<BuiltinSignedIntegral N> Int64(N n) : _n(n) { assert(_n==n); }
    int64_t get_si() const { return _n; }
};

uint32_t fac(uint8_t n);
uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint32_t bin(uint8_t n, uint8_t k);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);

struct ExactTag;
class Integer;
template<> struct IsNumber<Nat> : True { };
template<> struct IsNumber<Int> : True { };
template<> struct IsNumber<Integer> : True { };

//! \ingroup NumericModule
//! \brief Arbitrarily-sized integers.
class Integer
    : DeclareRingOperations<Integer,Integer,Natural>
    , DeclareLatticeOperations<Integer,Natural>
    , DeclareComparisonOperations<Integer,Boolean,Boolean>
    , DefineRingOperators<Integer>
    , DefineComparisonOperators<Integer,Boolean,Boolean>
{
  public:
    mpz_t _mpz;
  public:
    typedef ExactTag Paradigm;
  public:
    ~Integer();
    Integer();
    template<BuiltinUnsignedIntegral M> Integer(M m);
    template<BuiltinSignedIntegral N> Integer(N n);
    explicit Integer(const mpz_t);
    Integer(const Integer&);
    Integer(Integer&&);
    Integer& operator=(const Integer&);
    Integer& operator=(Integer&&);
    operator ExactNumber () const;

    friend Rational operator/(Integer const& z1, Integer const& z2);
    friend Integer& operator++(Integer& z);
    friend Integer& operator--(Integer& z);
    friend Integer& operator+=(Integer& z1, Integer const& z2);
    friend Integer& operator*=(Integer& z1, Integer const& z2);

    friend Int log2floor(Natural const& n);

    friend Dyadic hlf(Integer const& z);
    friend Rational rec(Integer const& z);
    friend Rational div(Integer const& z1, Integer const& z2);
    friend Rational pow(Integer const& z, Int n);
    friend Integer quot(Integer const& z1, Integer const& z2);
    friend Integer rem(Integer const& z1, Integer const& z2);
    friend Integer operator%(Integer const& z1, Integer const& z2);

    friend Bool is_nan(Integer const& z);
    friend Bool is_inf(Integer const& z);
    friend Bool is_finite(Integer const& z);
    friend Bool is_zero(Integer const& z);

    friend Sign sgn(Integer const& z);
    friend Comparison cmp(Integer const& z1, Integer const& z2);

    friend OutputStream& operator<<(OutputStream& os, Integer const& z);
    friend Integer operator"" _z(unsigned long long int n);
/*
    // Comparisons with arbitary ints go through Int64
    template<BuiltinIntegral N> friend inline auto operator==(Integer const& x, N n) -> decltype(x==Int64(n)) { return x==Int64(n); }
    template<BuiltinIntegral N> friend inline auto operator!=(Integer const& x, N n) -> decltype(x!=Int64(n)) { return x!=Int64(n); }
    template<BuiltinIntegral N> friend inline auto operator< (Integer const& x, N n) -> decltype(x!=Int64(n)) { return x< Int64(n); }
    template<BuiltinIntegral N> friend inline auto operator> (Integer const& x, N n) -> decltype(x!=Int64(n)) { return x> Int64(n); }
    template<BuiltinIntegral N> friend inline auto operator<=(Integer const& x, N n) -> decltype(x!=Int64(n)) { return x<=Int64(n); }
    template<BuiltinIntegral N> friend inline auto operator>=(Integer const& x, N n) -> decltype(x!=Int64(n)) { return x>=Int64(n); }
*/
  public:
    template<BuiltinIntegral N> N get() const;
    long int get_si() const;
    mpz_t const& get_mpz() const;
  private:
  public:
    Integer(Nat32 m);
    Integer(Int32 n);
    Integer(Nat64 m);
    Integer(Int64 n);
};

template<BuiltinUnsignedIntegral M> inline Integer::Integer(M m) : Integer(Nat64(m)) { }
template<BuiltinSignedIntegral N> inline Integer::Integer(N n) : Integer(Int64(n)) { }
Integer operator"" _z(unsigned long long int n);

template<BuiltinIntegral N> inline N Integer::get() const {
    N n=static_cast<N>(this->get_si()); ARIADNE_ASSERT(Integer(n)==*this); return n; }

template<class R, class A> R integer_cast(const A& a) {
    return cast_integer(a).template get<R>(); }

template<> class Positive<Integer> : public Integer {
  public:
    Positive<Integer>() : Integer() { }
    template<BuiltinUnsignedIntegral M> Positive<Integer>(M m) : Integer(m) { }
    Positive<Integer>(int n) = delete;
    explicit Positive<Integer>(Integer const& z) : Integer(z) { assert(z>=0); }
};

//! \brief A positive integer.
class Natural : public Positive<Integer> {
  public:
    Natural() : Positive<Integer>() { }
    template<BuiltinUnsignedIntegral M> Natural(M m) : Positive<Integer>(m) { }
    Natural(int n) = delete;
    explicit Natural(Integer const& z) : Positive<Integer>(z) { assert(z>=Integer(0)); }
    friend Natural& operator++(Natural& n) { ++static_cast<Integer&>(n); return n; }
    friend Natural& operator+=(Natural& n1, Natural const& n2) { static_cast<Integer&>(n1)+=n2; return n1; }
    friend Natural operator+(Natural const& n1, Natural const& n2) { return Natural(static_cast<Integer const&>(n1)+static_cast<Integer const&>(n2)); }
    friend Natural operator*(Natural const& n1, Natural const& n2) { return Natural(static_cast<Integer const&>(n1)*static_cast<Integer const&>(n2)); }
};

inline Natural cast_positive(Integer const& z) { return Natural(z); }

} // namespace Ariadne

#endif
