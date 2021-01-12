/***************************************************************************
 *            numeric/float_literals.hpp
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

/*! \file numeric/float_literals.hpp
 *  \brief Inclusion header for floating-point extended literals.
 */

#ifndef ARIADNE_FLOAT_LITERALS_HPP
#define ARIADNE_FLOAT_LITERALS_HPP

#include "float.decl.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Indicate that a floating-point literal only represents an double-precision approximation to a number.
ApproximateDouble operator"" _a (long double lx);
//! \ingroup NumericModule
//! \brief Indicate that a floating-point literal is the exact value of a number in double-precision.
//! \details The input \a lx is converted first converted to a double-precision value \a x.
//! If \a lx and \a x differ, then the decimal literal input to \a lx almost certainly does not represent a double-precision number exactly, and an assertion fails. <p/>
//! For example, <c>0.625_x</c> yields the exact double-precision value \f$5/2^3\f$, but <c>0.6_x</c> fails, since \f$3/5\f$ is not exactly representable as a double-precision number.
ExactDouble operator"" _x(long double lx);

FloatDPError operator"" _error(long double lx);
FloatDPBall operator"" _near(long double lx);
FloatDPUpperBound operator"" _upper(long double lx);
FloatDPLowerBound operator"" _lower(long double lx);
FloatDPApproximation operator"" _approx(long double lx);

} // namespace Ariadne

#endif
