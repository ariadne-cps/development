/***************************************************************************
 *            real.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file real.h
 *  \brief General (computable) real number class.
 */

#ifndef ARIADNE_REAL_H
#define ARIADNE_REAL_H

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <memory>

#include "utility/declarations.h"

#include "utility/tribool.h"
#include "numeric/rounding.h"
#include "utility/macros.h"
#include "numeric/integer.h"
#include "numeric/rational.h"

#include "numeric/float-approximate.h"
#include "numeric/float-validated.h"
#include "numeric/float-exact.h"

#include "geometry/interval.h"

namespace Ariadne {

class RealInterface;

//! \ingroup NumericModule
//! \brief Computable real numbers.
//!
//! Support over-approximation by an ExactInterval and approximation by a Float.
class Real {
    std::shared_ptr<RealInterface> _ptr;
  public:
    typedef Real NumericType;
  public:
    explicit Real(RealInterface* raw_ptr);
  public:
    //! Destructor.
    ~Real();
    //! Default constructor yields the exact value \a 0.
    Real();
    //! \brief Construct from a string literal.
    //! This can be used to create real numbers representing decimal values e.g. \c %Real("4.2") .
    explicit Real(const std::string& s);
    //! \brief Construct from double-precision values giving lower and upper bounds for the exact value.
    explicit Real(double l, double u);
    //! \brief Construct from double-precision values giving lower and upper bounds for the exact value and a nearest approximation.
    explicit Real(double l, double x, double u);

    //! \brief Convert from a builtin natural number.
    Real(unsigned int m);
    //! \brief Convert from a builtin integer.
    Real(int n);
    //! \brief Convert from a builtin double-precision floating-point value.
    //! A numeric literal is first processed by the language support, and the resulting %Real may not have
    //! the same value as the mathematical literal. e.g. \c %Real(4.2) has the value 4.2000000000000002 to 16 decimal places.
    explicit Real(double x);
#ifdef HAVE_GMPXX_H
    //! \brief Convert from an integer.
    Real(const Integer& z);
    //! \brief Convert from a rational number.
    Real(const Rational& q);
#endif
    //! \brief Convert from a decimal number.
    Real(const Decimal& d);
    //! \brief Convert from a dyadic number.
    Real(const Dyadic& x);
    //! \brief Convert from an exact floating-point number.
    Real(const ExactFloat& x);
    //! \brief Copy constructor.
    Real(const Real&);
    //! \brief Copy assignment.
    Real& operator=(const Real&);

    explicit operator ApproximateFloat() const;
    explicit operator ValidatedFloat() const;
    operator UpperFloat() const;
  public:
    //! \brief Get an approximation as a builtin double-precision floating-point number.
    double get_d() const;
  private:
    friend class ApproximateFloat;
    friend class ValidatedFloat;
    friend OutputStream& operator<<(OutputStream& os, const Real& x);
};

//@{
//! \related Real \name Arithmetic Operators

//!  \brief Unary plus operator.
Real operator+(const Real& x);
//!  \brief Unary negation operator.
Real operator-(const Real& x);
//!  \brief Binary addition operator.
Real operator+(const Real& x, const Real& y);
//!  \brief Subtraction operator.
Real operator-(const Real& x, const Real& y);
//! \brief Binary multiplication operator.
Real operator*(const Real& x, const Real& y);
//! \brief Division operator.
Real operator/(const Real& x, const Real& y);

//!  \brief Inplace addition operator.
inline Real& operator+=(Real& x, const Real& y) { return x=x+y; }
//!  \brief Inplace subtraction operator.
inline Real& operator-=(Real& x, const Real& y) { return x=x-y; }
//!  \brief Inplace multiplication operator.
inline Real& operator*=(Real& x, const Real& y) { return x=x*y; }
//!  \brief Inplace division operator.
inline Real& operator/=(Real& x, const Real& y) { return x=x/y; }

//@}

template<class N, typename std::enable_if<std::is_integral<N>::value,int>::type=0>
    inline Real operator*(const Real& i1, N n2) { return operator*(i1,Real(n2)); };
template<class N, typename std::enable_if<std::is_integral<N>::value,int>::type=0>
    inline Real operator*(N n1, const Real& i2) { return operator*(Real(n1),i2); };
template<class N, typename std::enable_if<std::is_integral<N>::value,int>::type=0>
    inline Real operator/(const Real& i1, N n2) { return operator/(i1,Real(n2)); };
template<class N, typename std::enable_if<std::is_integral<N>::value,int>::type=0>
    inline Real operator/(N n1, const Real& i2) { return operator/(Real(n1),i2); };

/*
inline Real operator+(const Real& x, double y) { return x+static_cast<Real>(y); }
inline Real operator-(const Real& x, double y) { return x-static_cast<Real>(y); }
inline Real operator*(const Real& x, double y) { return x*static_cast<Real>(y); }
inline Real operator/(const Real& x, double y) { return x/static_cast<Real>(y); }
inline Real operator+(double x, const Real& y) { return static_cast<Real>(x)+y; }
inline Real operator-(double x, const Real& y) { return static_cast<Real>(x)-y; }
inline Real operator*(double x, const Real& y) { return static_cast<Real>(x)*y; }
inline Real operator/(double x, const Real& y) { return static_cast<Real>(x)/y; }
*/

/*
inline ApproximateFloat operator+(const Real& x, const ApproximateFloat& y) { return static_cast<ApproximateFloat>(x)+y; }
inline ApproximateFloat operator-(const Real& x, const ApproximateFloat& y) { return static_cast<ApproximateFloat>(x)-y; }
inline ApproximateFloat operator*(const Real& x, const ApproximateFloat& y) { return static_cast<ApproximateFloat>(x)*y; }
inline ApproximateFloat operator/(const Real& x, const ApproximateFloat& y) { return static_cast<ApproximateFloat>(x)/y; }
inline ApproximateFloat operator+(const ApproximateFloat& x, const Real& y) { return x+static_cast<ApproximateFloat>(y); }
inline ApproximateFloat operator-(const ApproximateFloat& x, const Real& y) { return x-static_cast<ApproximateFloat>(y); }
inline ApproximateFloat operator*(const ApproximateFloat& x, const Real& y) { return x*static_cast<ApproximateFloat>(y); }
inline ApproximateFloat operator/(const ApproximateFloat& x, const Real& y) { return x/static_cast<ApproximateFloat>(y); }
inline ApproximateFloat& operator+=(ApproximateFloat& x, const Real& y) { return x+=static_cast<ApproximateFloat>(y); }
inline ApproximateFloat& operator-=(ApproximateFloat& x, const Real& y) { return x-=static_cast<ApproximateFloat>(y); }
inline ApproximateFloat& operator*=(ApproximateFloat& x, const Real& y) { return x*=static_cast<ApproximateFloat>(y); }
inline ApproximateFloat& operator/=(ApproximateFloat& x, const Real& y) { return x/=static_cast<ApproximateFloat>(y); }

inline ValidatedFloat operator+(const Real& x, const ValidatedFloat& y) { return static_cast<ValidatedFloat>(x)+y; }
inline ValidatedFloat operator-(const Real& x, const ValidatedFloat& y) { return static_cast<ValidatedFloat>(x)-y; }
inline ValidatedFloat operator*(const Real& x, const ValidatedFloat& y) { return static_cast<ValidatedFloat>(x)*y; }
inline ValidatedFloat operator/(const Real& x, const ValidatedFloat& y) { return static_cast<ValidatedFloat>(x)/y; }
inline ValidatedFloat operator+(const ValidatedFloat& x, const Real& y) { return x+static_cast<ValidatedFloat>(y); }
inline ValidatedFloat operator-(const ValidatedFloat& x, const Real& y) { return x-static_cast<ValidatedFloat>(y); }
inline ValidatedFloat operator*(const ValidatedFloat& x, const Real& y) { return x*static_cast<ValidatedFloat>(y); }
inline ValidatedFloat operator/(const ValidatedFloat& x, const Real& y) { return x/static_cast<ValidatedFloat>(y); }
inline ValidatedFloat& operator+=(ValidatedFloat& x, const Real& y) { return x+=static_cast<ValidatedFloat>(y); }
inline ValidatedFloat& operator-=(ValidatedFloat& x, const Real& y) { return x-=static_cast<ValidatedFloat>(y); }
inline ValidatedFloat& operator*=(ValidatedFloat& x, const Real& y) { return x*=static_cast<ValidatedFloat>(y); }
inline ValidatedFloat& operator/=(ValidatedFloat& x, const Real& y) { return x/=static_cast<ValidatedFloat>(y); }

inline ExactInterval operator+(const Real& x, const ExactInterval& y) { return static_cast<ExactInterval>(x)+y; }
inline ExactInterval operator-(const Real& x, const ExactInterval& y) { return static_cast<ExactInterval>(x)-y; }
inline ExactInterval operator*(const Real& x, const ExactInterval& y) { return static_cast<ExactInterval>(x)*y; }
inline ExactInterval operator/(const Real& x, const ExactInterval& y) { return static_cast<ExactInterval>(x)/y; }
inline ExactInterval operator+(const ExactInterval& x, const Real& y) { return x+static_cast<ExactInterval>(y); }
inline ExactInterval operator-(const ExactInterval& x, const Real& y) { return x-static_cast<ExactInterval>(y); }
inline ExactInterval operator*(const ExactInterval& x, const Real& y) { return x*static_cast<ExactInterval>(y); }
inline ExactInterval operator/(const ExactInterval& x, const Real& y) { return x/static_cast<ExactInterval>(y); }
inline ExactInterval& operator+=(ExactInterval& x, const Real& y) { return x+=static_cast<ExactInterval>(y); }
inline ExactInterval& operator-=(ExactInterval& x, const Real& y) { return x-=static_cast<ExactInterval>(y); }
inline ExactInterval& operator*=(ExactInterval& x, const Real& y) { return x*=static_cast<ExactInterval>(y); }
inline ExactInterval& operator/=(ExactInterval& x, const Real& y) { return x/=static_cast<ExactInterval>(y); }
*/

//@{
//! \related Real \name Comparison Operators

//! \brief Equality operator.
//! Returns \c true if the numbers have the same representation or can be evaluated exactly to the same value.
//! Returns \c false if the numbers can be evalated to different values.
//! Returns \c indeterminate if the numbers cannot be shown to be the same or different. Implementation dependent.
inline tribool operator==(const Real& x, const Real& y) { return static_cast<ExactInterval>(x)==static_cast<ExactInterval>(y); }
//! \brief Inequality operator.
//! Returns \c true if the numbers have been proved to be different, such as by evaluation to different values.
//! Returns \c false if the numbers have the same representation or can be evaluated exactly to the same value.
//! Returns \c indeterminate if the numbers cannot be shown to be the same or different. Implementation dependent.
inline tribool operator!=(const Real& x, const Real& y) { return static_cast<ExactInterval>(x)!=static_cast<ExactInterval>(y); }
//! \brief Greater-than-or-equal-to comparison operator.
inline tribool operator>=(const Real& x, const Real& y) { return static_cast<ExactInterval>(x)>=static_cast<ExactInterval>(y); }
//! \brief Less-than-or-equal-to comparison operator.
inline tribool operator<=(const Real& x, const Real& y) { return static_cast<ExactInterval>(x)<=static_cast<ExactInterval>(y); }
//! \brief Strictly-greater-than comparison operator.
inline tribool operator> (const Real& x, const Real& y) { return static_cast<ExactInterval>(x)> static_cast<ExactInterval>(y); }
//! \brief Strictly-less-than comparison operator.
inline tribool operator< (const Real& x, const Real& y) { return static_cast<ExactInterval>(x)< static_cast<ExactInterval>(y); }
//@}

inline tribool operator==(const Real& x, double y) { return static_cast<ExactInterval>(x)==static_cast<ExactInterval>(y); }
inline tribool operator!=(const Real& x, double y) { return static_cast<ExactInterval>(x)!=static_cast<ExactInterval>(y); }
inline tribool operator>=(const Real& x, double y) { return static_cast<ExactInterval>(x)>=static_cast<ExactInterval>(y); }
inline tribool operator<=(const Real& x, double y) { return static_cast<ExactInterval>(x)<=static_cast<ExactInterval>(y); }
inline tribool operator> (const Real& x, double y) { return static_cast<ExactInterval>(x)> static_cast<ExactInterval>(y); }
inline tribool operator< (const Real& x, double y) { return static_cast<ExactInterval>(x)< static_cast<ExactInterval>(y); }

//@{
//! \related Real \name Arithmetical, algebraic and transcendental functions

//!  \brief A floating-point bound for |x|.
Float mag(const Real& x);

//!  \brief The constant infinity.
extern const Real infinity;
//!  \brief The constant pi.
extern const Real pi;

//!  \brief The absolute value function \c |x|.
Real abs(const Real& x);
//!  \brief The unary plus function \c +x.
Real pos(const Real& x);
//!  \brief The unary negation function \c -x.
Real neg(const Real& x);
//!  \brief The reciprocal function \c 1/x.
Real sqr(const Real& x);
//!  \brief The reciprocal function \c 1/x.
Real rec(const Real& x);
//!  \brief The binary addition function \c x+y.
Real add(const Real& x, const Real& y);
//!  \brief The subtraction function \c x-y.
Real sub(const Real& x, const Real& y);
//!  \brief The binary multiplication function \c x*y.
Real mul(const Real& x, const Real& y);
//!  \brief The division function \c x/y.
Real div(const Real& x, const Real& y);
//!  \brief The positive integer power function \c x^m.
Real pow(const Real& x, uint m);
//!  \brief The integer power function \c x^n.
Real pow(const Real& x, int n);
//!  \brief The square-root function.
Real sqrt(const Real& x);
//!  \brief The exponential function.
Real exp(const Real& x);
//!  \brief The natural logarithm function.
Real log(const Real& x);
//!  \brief The sine function.
Real sin(const Real& x);
//!  \brief The cosine function.
Real cos(const Real& x);
//!  \brief The tangent function.
Real tan(const Real& x);
//!  \brief The inverse sine function.
Real asin(const Real& x);
//!  \brief The inverse cosine function.
Real acos(const Real& x);
//!  \brief The inverse tangent function.
Real atan(const Real& x);

//@}

//@{
//! \related Real \name Input/output operators

//! \related Real \brief Write to an output stream.
std::ostream& operator<<(std::ostream& os, const Real& x);
//@}


} // namespace Ariadne

#endif