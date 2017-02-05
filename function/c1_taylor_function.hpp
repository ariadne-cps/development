/***************************************************************************
 *            c1_taylor_function.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file c1_taylor_function.hpp
 *  \brief Over-approximations of continuously-differentiable functions based on Taylor expansions.
 */

#ifndef ARIADNE_C1_TAYLOR_FUNCTION_HPP
#define ARIADNE_C1_TAYLOR_FUNCTION_HPP

#include <iosfwd>
#include "utility/declarations.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "numeric/float.hpp"
#include "algebra/expansion.hpp"
#include "function/domain.hpp"

namespace Ariadne {


class MultiIndex;
class C1TaylorSeries;
class C1TaylorFunction;


/*! \ingroup FunctionModelSubModule
 *  \brief A C1TaylorSeries is a univariate function \f$f:\R\rightarrow\R\f$ on an interval \f$[a,b]\f$ is approximated by polynomial \f$p\f$ with error bounds \f$e_0 \geq |f(c)-p(c)|\f$ and \f$e_1\geq sup_{\xi\in D} |f'(\xi)-p'(\xi)|\f$.
 */
class C1TaylorSeries
{
  public:
    std::vector<Float64> _coefficients;
    Float64 _zero_error;
    Float64 _uniform_error;
    Float64 _derivative_error;
  private:
    C1TaylorSeries(Nat d);
  public:
    C1TaylorSeries();
    static C1TaylorSeries constant(Float64);
    static C1TaylorSeries coordinate();
  public:
    IntervalDomainType domain() const;
    Nat degree() const;
    Void sweep(Float64 threshold);
  public:
    friend C1TaylorSeries& operator+=(C1TaylorSeries&, ValidatedNumericType);
    friend C1TaylorSeries& operator*=(C1TaylorSeries&, ValidatedNumericType);
    friend C1TaylorSeries operator+(C1TaylorSeries, C1TaylorSeries);
    friend C1TaylorSeries operator*(C1TaylorSeries, C1TaylorSeries);
    friend ValidatedNumericType evaluate(C1TaylorSeries, ValidatedNumericType);
    friend C1TaylorSeries compose(C1TaylorSeries, C1TaylorSeries);
    friend OutputStream& operator<< (OutputStream& os, const C1TaylorSeries& f);
};

class C1TaylorFunction
{
  public:
    typedef Float64Bounds NumericType;
  public:
    Expansion<Float64> _expansion;
    Float64 _zero_error;
    Float64 _uniform_error;
    Array<Float64> _derivative_errors;
  private:
  public:
    C1TaylorFunction();
    C1TaylorFunction(SizeType as);
  public:
    static C1TaylorFunction constant(SizeType as, Float64 c);
    static C1TaylorFunction coordinate(SizeType as, SizeType ind);
  public:
    BoxDomainType domain() const;
    Nat argument_size() const;
    Void sweep(Float64 threshold);
    C1TaylorFunction& operator=(NumericType c);
    Void clear();
  public:
    friend C1TaylorFunction& operator+=(C1TaylorFunction& f, Float64 c);
    friend C1TaylorFunction& operator*=(C1TaylorFunction& f, Float64 c);
    friend C1TaylorFunction operator+(C1TaylorFunction f1, C1TaylorFunction f2);
    friend C1TaylorFunction operator*(C1TaylorFunction f1, C1TaylorFunction f2);
    friend NumericType evaluate(C1TaylorFunction f, Vector<NumericType> x);
    friend C1TaylorFunction compose(C1TaylorSeries f, C1TaylorFunction g);
    friend C1TaylorFunction compose(C1TaylorFunction f, Vector<C1TaylorFunction> g);
    friend OutputStream& operator<< (OutputStream& os, const C1TaylorFunction& f);
};

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_HPP
