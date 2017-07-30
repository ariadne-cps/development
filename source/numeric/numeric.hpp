/***************************************************************************
 *            numeric.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file numeric.hpp
 *  \brief Number classes. File suitable for use as a pre-compiled header.
 */

#ifndef ARIADNE_NUMERIC_HPP
#define ARIADNE_NUMERIC_HPP

#include "utility/standard.hpp"

#include "config.h"

#include "utility/declarations.hpp"

#include "numeric/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/decimal.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float.hpp"
#include "numeric/real.hpp"
#include "numeric/number.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline Int numeric_cast(const Float64& a) { return Int(a.get_d()); }
template<> inline Int numeric_cast(const FloatMP& a) { return Int(a.get_d()); }
template<> inline double numeric_cast(const Float64& a) { return a.get_d(); }
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline double numeric_cast(const Float64Value& a) { return a.get_d(); }
template<> inline double numeric_cast(const Float64Bounds& a) { return a.get_d(); }
template<> inline double numeric_cast(const Float64Approximation& a) { return a.get_d(); }
template<> inline float numeric_cast(const double& a) { return a; }
template<> inline float numeric_cast(const Float64& a) { return a.get_d(); }
template<> inline float numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Float64 numeric_cast(const Float64Value& a) { return a.raw(); }

template<> inline Real numeric_cast(const Float64& a) { return Real(ExactDouble(a.get_d())); }
template<> inline Real numeric_cast(const Float64Value& a) { return numeric_cast<Real>(a.raw()); }
template<> inline Real numeric_cast(const Float64Bounds& a) { return numeric_cast<Real>(Float64Approximation(a).raw()); }

template<> inline Float64Ball numeric_cast(const Real& a) { return Float64Ball(a,Precision64()); }
template<> inline Float64Bounds numeric_cast(const Real& a) { return Float64Bounds(a,Precision64()); }
template<> inline Float64Approximation numeric_cast(const Real& a) { return Float64Approximation(a,Precision64()); }

} // namespace Ariadne

#endif