/***************************************************************************
 *            numeric/integer.cpp
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

/*! \file numeric/integer.cpp
 *  \brief
 */


#include "utility/stdlib.hpp"

#include "integer.hpp"

#include "utility/macros.hpp"
#include "utility/string.hpp"
#include "numeric/logical.hpp"

#include <limits>

namespace Ariadne {

Comparison cmp(Integer const& z1, Integer const& z2);
Integer make_integer(unsigned long long int n);

uint32_t
fac(uint8_t n)
{
    ARIADNE_ASSERT(n<13); // Maximum factorial in 32 bits
    uint32_t  r=1;
    for(uint8_t i=1; i<=n; ++i) {
        r*=i;
    }
    return r;
}

uint32_t
bin(uint8_t n, uint8_t k)
{
    ARIADNE_ASSERT(n<32);  // Maximum computable bin(n,n/2) using 32 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { ARIADNE_FAIL_MSG("bin("<<n<<","<<k<<")"); }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint32_t r=1;
    for(uint8_t i=1; i<=k; ++i) {
        r*=(n+1u-i);
        r/=i;
    }
    return r;
}

uint16_t
fac(uint16_t n)
{
    ARIADNE_ASSERT(n<9); // Maximum factorial in 16 bits
    uint16_t  r=1;
    for(uint16_t i=1; i<=n; ++i) {
        r*=i;
    }
    return r;
}


uint16_t
bin(uint16_t n, uint16_t k)
{
    ARIADNE_ASSERT(n<16);  // Maximum computable bin(n,n/2) using 16 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { ARIADNE_FAIL_MSG("bin("<<n<<","<<k<<")"); }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint16_t r=1;
    for(uint16_t i=1; i<=k; ++i) {
        r*=(n+1-i);
        r/=i;
    }
    return r;
}

uint32_t
fac(uint32_t n)
{
    ARIADNE_ASSERT(n<13); // Maximum factorial in 32 bits
    uint32_t  r=1;
    for(uint32_t i=1; i<=n; ++i) {
        r*=i;
    }
    return r;
}


uint32_t
bin(uint32_t n, uint32_t k)
{
    ARIADNE_ASSERT(n<31);  // Maximum computable bin(n,n/2) using 32 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { ARIADNE_FAIL_MSG("bin("<<n<<","<<k<<")"); }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint32_t r=1;
    for(uint32_t i=1; i<=k; ++i) {
        r*=(n+1-i);
        r/=i;
    }
    return r;
}



uint64_t
fac(uint64_t n)
{
    ARIADNE_ASSERT(n<21); // Maximum factorial in 64 bits
    uint64_t  r=1;
    for(uint64_t i=1; i<=n; ++i) {
        r*=i;
    }
    return r;
}


uint64_t
bin(uint64_t n, uint64_t k)
{
    ARIADNE_ASSERT(n<63);  // Maximum computable bin(n,n/2) using 64 bits
                           // Note that this is shorter than the maximum representable factorial
    if(k>n+1) { ARIADNE_FAIL_MSG("bin("<<n<<","<<k<<")"); }
    if(k==n+1) { return 0; }
    ARIADNE_ASSERT(k<=n);
    uint64_t r=1;
    for(uint64_t i=1; i<=k; ++i) {
        r*=(n+1-i);
        r/=i;
    }
    return r;
}


Integer::~Integer() {
    mpz_clear(_mpz);
}

Integer::Integer() {
    mpz_init(_mpz);
    mpz_set_si(_mpz,0);
}

Integer::Integer(Nat32 m) {
    mpz_init(_mpz);
    mpz_set_ui(_mpz,m.get_ui());
}

Integer::Integer(Int32 n) {
    mpz_init(_mpz);
    mpz_set_ui(_mpz,static_cast<long unsigned int>(n.get_si()));
}

Integer::Integer(Nat64 m) {
    mpz_init(_mpz);
    static const unsigned int max=std::numeric_limits<Int>::max();
    static const unsigned long long int lmax=max;
    static const Integer zmax=Integer(Int32(max));
    unsigned long long int larg = m.get_ui();
    unsigned long long int lquot = larg / lmax;
    unsigned long long int lrem = larg % lmax;
    unsigned int rem = static_cast<unsigned int>(lrem);
    unsigned int quot = static_cast<unsigned int>(lquot);
    assert(larg==lquot*lmax+lrem);
    mpz_set_ui(_mpz,rem);
    if(lquot!=0) {
        *this += Integer(quot)*zmax;
    }
}

Integer::Integer(Int64 larg) {
    mpz_init(_mpz);
    if(larg.get_si()<0) {
        *this = -Integer(Nat64(-larg.get_si()));
    } else {
        *this = Integer(Nat64(larg.get_si()));
    }
}

Integer::Integer(const mpz_t z) {
    mpz_init(_mpz);
    mpz_set(_mpz,z);
}

Integer::Integer(const Integer& z) {
    mpz_init(_mpz);
    mpz_set(_mpz,z._mpz);
}

Integer::Integer(Integer&& z) {
    mpz_init(_mpz);
    mpz_swap(_mpz,z._mpz);
}

Integer& Integer::operator=(const Integer& z) {
    mpz_set(_mpz,z._mpz);
    return *this;
}

Integer& Integer::operator=(Integer&& z) {
    mpz_swap(_mpz,z._mpz);
    return *this;
}


mpz_t const& Integer::get_mpz() const {
    return this->_mpz;
}

long int Integer::get_si() const {
    return mpz_get_si(this->_mpz);
}


Integer& operator++(Integer& z) {
    mpz_add_ui(z._mpz,z._mpz,1u);
    return z;
}

Integer& operator--(Integer& z) {
    mpz_sub_ui(z._mpz,z._mpz,1u);
    return z;
}

/*
Integer& operator+=(Integer& z1, Integer const& z2) {
    mpz_add(z1._mpz,z1._mpz,z2._mpz);
    return z1;
}

Integer& operator-=(Integer& z1, Integer const& z2) {
    mpz_sub(z1._mpz,z1._mpz,z2._mpz);
    return z1;
}

Integer& operator*=(Integer& z1, Integer const& z2) {
    mpz_mul(z1._mpz,z1._mpz,z2._mpz);
    return z1;
}
*/

Integer nul(Integer const& z) {
    Integer r;
    mpz_set_si(r._mpz,0);
    return r;
}

Integer pos(Integer const& z) {
    Integer r;
    mpz_set(r._mpz,z._mpz);
    return r;
}

Integer neg(Integer const& z) {
    Integer r;
    mpz_neg(r._mpz,z._mpz);
    return r;
}

Natural sqr(Integer const& z) {
    Natural r;
    mpz_mul(r._mpz,z._mpz,z._mpz);
    return r;
}

Integer add(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_add(r._mpz,z1._mpz,z2._mpz);
    return r;
}

Integer sub(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_sub(r._mpz,z1._mpz,z2._mpz);
    return r;
}

Integer mul(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_mul(r._mpz,z1._mpz,z2._mpz);
    return r;
}

Integer quot(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_div(r._mpz,z1._mpz,z2._mpz);
    return r;
}

Integer rem(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_cdiv_r(r._mpz,z1._mpz,z2._mpz);
    return r;
}

Integer pow(Integer const& z, Nat m) {
    unsigned long int lm=m;
    Integer r;
    mpz_pow_ui(r._mpz,z._mpz,lm);
    return r;
}


Integer min(Integer const& z1,Integer const& z2) {
    return (z1<z2)?z1:z2;
}

Integer max(Integer const& z1,Integer const& z2) {
    return (z1>z2)?z1:z2;
}

Natural abs(Integer const& z) {
    Natural r;
    mpz_abs(r._mpz,z._mpz);
    return r;
}

Natural max(Natural const& z1,Natural const& z2) {
    return (z1>z2)?z1:z2;
}

Natural min(Natural const& z1,Natural const& z2) {
    return (z1<z2)?z1:z2;
}


Bool is_nan(Integer const& z) {
    return false;
}

Bool is_inf(Integer const& z) {
    return false;
}

Bool is_finite(Integer const& z) {
    return true;
}

Bool is_zero(Integer const& z) {
    return mpz_cmp_si(z._mpz,0)==0;
}

Sign sgn(Integer const& z) {
    return static_cast<Sign>(mpz_sgn(z._mpz));
}

Comparison cmp(Integer const& z1, Integer const& z2) {
    auto c=mpz_cmp(z1._mpz,z2._mpz);
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

Boolean eq(Integer const& z1, Integer const& z2) {
    return mpz_cmp(z1._mpz,z2._mpz)==0;
}

Boolean lt(Integer const& z1, Integer const& z2) {
    return mpz_cmp(z1._mpz,z2._mpz) < 0;
}


//   mpz_get_str (char *str, mpz_exp_t *expptr, Int b, SizeType n, mpz_t op, mpz_rnd_t rnd)
// If str is not a null pointer, it should point to a block of storage large enough for the significand,
// i.e., at least maq1(n + 2, 7). The extra two bytes are for a possible minus sign,
// and for the terminating null character, and the value 7 accounts for -@Inf@ plus the terminating null character.
OutputStream& operator<<(OutputStream& os, Integer const& z) {
    static const int OUTPUT_BUFFER_SIZE=4096;
    char str[OUTPUT_BUFFER_SIZE];
    str[OUTPUT_BUFFER_SIZE-1]='\0';
    mpz_get_str (str, 10, z._mpz);
    assert(str[OUTPUT_BUFFER_SIZE-1]=='\0');
    return os << str;
}

Integer make_integer(unsigned long long int n) {
    static const unsigned int max=std::numeric_limits<Int>::max();
    static const unsigned long long int m=max;
    unsigned long long int q = n / m;
    unsigned long long int r = n % m;
    unsigned int rem = static_cast<unsigned int>(r);
    assert(n==q*m+r);
    if(q==0) {
        return Integer(rem);
    } else {
        return make_integer(q)*Integer(max)+Integer(rem);
    }
}

Integer operator"" _z(unsigned long long int n) {
    return Integer(Nat64(n));
}


template<> inline String class_name<uint>() { return "uint"; }

template<> inline String class_name<int>() { return "int"; }

template<> String class_name<Integer>() { return "Integer"; }

template<> String class_name<Natural>() { return "Natural"; }

Int log2floor(Natural const& z) {
    return mpz_sizeinbase(z._mpz,2)-1;
}

OutputStream& operator<<(OutputStream& os, Sign s) {
    return os << ( (s==Sign::ZERO) ? "ZERO" : (s==Sign::NEGATIVE) ? "NEGATIVE" : "POSITIVE" );
}

static_assert(IsOrderedLatticeRing<Integer>);
static_assert(AreMixedOrderedLatticeRing<Integer,Int>);


} // namespace Ariadne
