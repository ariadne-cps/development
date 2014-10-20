/***************************************************************************
 *            integer.cc
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file integer.cc
 *  \brief
 */



#include "utility/stdlib.h"

#include "integer.h"

#include "utility/string.h"
#include "numeric/logical.h"

#include <limits>

namespace Ariadne {

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
    mpz_set_ui(_mpz,n.get_si());
}

Integer::Integer(Nat64 m) {
    mpz_init(_mpz);
    static const unsigned int max=std::numeric_limits<int>::max();
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

Integer operator+(Integer const& z) {
    return Integer(z);
}

Integer operator-(Integer const& z) {
    Integer r; mpz_neg(r._mpz,z._mpz);
    return std::move(r);
}

Integer operator+(Integer const& z1, Integer const& z2) {
    Integer r; mpz_add(r._mpz,z1._mpz,z2._mpz);
    return std::move(r);
}

Integer operator-(Integer const& z1, Integer const& z2) {
    Integer r; mpz_sub(r._mpz,z1._mpz,z2._mpz);
    return std::move(r);
}

Integer operator*(Integer const& z1, Integer const& z2) {
    Integer r; mpz_mul(r._mpz,z1._mpz,z2._mpz);
    return std::move(r);
}

Integer& operator++(Integer& z) {
    mpz_add_ui(z._mpz,z._mpz,1u);
    return z;
}

Integer& operator--(Integer& z) {
    mpz_sub_ui(z._mpz,z._mpz,1u);
    return z;
}

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

Integer pos(Integer const& z) {
    Integer r; mpz_set(r._mpz,z._mpz);
    return std::move(r);
}

Integer neg(Integer const& z) {
    Integer r; mpz_neg(r._mpz,z._mpz);
    return std::move(r);
}

Integer sqr(Integer const& z) {
    return z*z;
}

Integer add(Integer const& z1, Integer const& z2) {
    Integer r; mpz_add(r._mpz,z1._mpz,z2._mpz);
    return std::move(r);
}

Integer sub(Integer const& z1, Integer const& z2) {
    Integer r; mpz_sub(r._mpz,z1._mpz,z2._mpz);
    return std::move(r);
}

Integer mul(Integer const& z1, Integer const& z2) {
    Integer r; mpz_mul(r._mpz,z1._mpz,z2._mpz);
    return std::move(r);
}

Integer pow(Integer const& z, Nat m) {
    unsigned long int lm=m;
    Integer r;
    mpz_pow_ui(r._mpz,z._mpz,lm);
    return std::move(r);
}


Integer abs(Integer const& z) {
    Integer r; mpz_abs(r._mpz,z._mpz);
    return std::move(r);
}

Integer const& min(Integer const& z1,Integer const& z2) {
    return (z1<z2)?z1:z2;
}

Integer const& max(Integer const& z1,Integer const& z2) {
    return (z1>z2)?z1:z2;
}

OutputStream& operator<<(OutputStream& os, Comparison const& cmp) {
    switch(cmp) {
        case Comparison::LESS: os << "LESS";  break;
        case Comparison::EQUAL: os << "EQUAL";  break;
        case Comparison::GREATER: os << "GREATER";  break;
    }
    return os;
}

Comparison cmp(Integer const& z1, Integer const& z2) {
    auto c=mpz_cmp(z1._mpz,z2._mpz);
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

Boolean operator==(Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)==Comparison::EQUAL;
}

Boolean operator!=(Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)!=Comparison::EQUAL;
}

Boolean operator<=(Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)!=Comparison::GREATER;
}

Boolean operator>=(Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)!=Comparison::LESS;
}

Boolean operator< (Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)==Comparison::LESS;
}

Boolean operator> (Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)==Comparison::GREATER;
}

//   mpz_get_str (char *str, mpz_exp_t *expptr, int b, size_t n, mpz_t op, mpz_rnd_t rnd)
// If str is not a null pointer, it should point to a block of storage large enough for the significand,
// i.e., at least maq1(n + 2, 7). The extra two bytes are for a possible minus sign,
// and for the terminating null character, and the value 7 accounts for -@Inf@ plus the terminating null character.
OutputStream& operator<<(OutputStream& os, Integer const& z) {
    char str[255];
    str[254]='\0';
    mpz_get_str (str, 10, z._mpz);
    assert(str[254]=='\0');
    return os << str;
}

Integer make_integer(unsigned long long int n) {
    static const unsigned int max=std::numeric_limits<int>::max();
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

template<> String class_name<Integer>() { return "Integer"; }

} // namespace Ariadne
