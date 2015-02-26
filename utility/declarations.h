/***************************************************************************
 *            declarations.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file declarations.h
 *  \brief Forward declarations of types and classes.
 */

#ifndef ARIADNE_DECLARATIONS_H
#define ARIADNE_DECLARATIONS_H

#include <iosfwd>

#include "utility/metaprogramming.h"
#include "utility/typedefs.h"

#include "numeric/paradigm.h"
#include "numeric/logical.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"

#include "geometry/interval.decl.h"
#include "geometry/box.decl.h"

namespace Ariadne {

typedef unsigned int uint;

//! Internal name for output stream.
typedef OutputStream OutputStream;

//! Internal name for void type.
typedef void Void;

//! Internal name for builtin boolean type.
typedef bool Bool;
//! Internal name for builtin integers.
typedef int Int;
//! Internal name for builtin unsigned integers.
typedef uint Nat;

// Define as a class for consistency with other value types
class String;
class Kleenean;

typedef SizeType DimensionType;

typedef ErrorFloat64 ValidatedNormType;
typedef ApproximateFloat64 ApproximateNormType;

typedef ErrorFloat64 NormType; // FIXME: Remove this typedef
typedef ErrorFloat64 ErrorType; // FIXME: Remove this typedef
typedef ApproximateFloat64 ApproximateErrorType; // FIXME: Remove this typedef

template<class I> struct CanonicalNumberTypedef;
template<> struct CanonicalNumberTypedef<ExactTag> { typedef ExactNumber Type; };
template<> struct CanonicalNumberTypedef<EffectiveTag> { typedef EffectiveNumber Type; };
template<> struct CanonicalNumberTypedef<ValidatedTag> { typedef ValidatedNumber Type; };
template<> struct CanonicalNumberTypedef<ApproximateTag> { typedef ApproximateNumber Type; };
template<class I> using CanonicalNumberType = typename CanonicalNumberTypedef<I>::Type;

template<class X> struct InformationTypedef;
template<> struct InformationTypedef<ExactNumber> { typedef ExactTag Type; };
template<> struct InformationTypedef<EffectiveNumber> { typedef EffectiveTag Type; };
template<> struct InformationTypedef<ValidatedNumber> { typedef ValidatedTag Type; };
template<> struct InformationTypedef<ApproximateNumber> { typedef ApproximateTag Type; };
template<class X> using InformationTag = typename InformationTypedef<X>::Type;

template<class X> using Scalar = X;
// Concrete class declarations
template<class X> class Vector;
template<class X> class Covector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class UnivariateDifferential;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class Series;

template<class X> class AffineModel;
template<class P, class F> class TaylorModel;
template<class X> class Formula;
template<class X> class Algebra;

template<class X> class Point;

typedef Vector<Rational> RationalVector;
typedef Vector<Real> RealVector;
typedef Vector<Float64> FloatVector;
typedef Vector<RawFloat64> RawFloatVector;
typedef Vector<ApproximateFloat64> ApproximateFloatVector;
typedef Vector<ValidatedFloat64> ValidatedFloatVector;
typedef Vector<ExactFloat64> ExactFloatVector;

typedef Vector<ApproximateNumber> ApproximateVector;
typedef Vector<ValidatedNumber> ValidatedVector;
typedef Vector<EffectiveNumber> EffectiveVector;
typedef Vector<ExactNumber> ExactVector;

typedef Matrix<Rational> RationalMatrix;
typedef Matrix<Real> RealMatrix;
typedef Matrix<RawFloat64> RawFloatMatrix;
typedef Matrix<Float64> FloatMatrix;
typedef Matrix<ApproximateFloat64> ApproximateFloatMatrix;
typedef Matrix<ValidatedFloat64> ValidatedFloatMatrix;
typedef Matrix<ExactFloat64> ExactFloatMatrix;

typedef Point<ApproximateNumber> ApproximatePoint;
typedef Point<ValidatedNumber> ValidatedPoint;
typedef Point<EffectiveNumber> EffectivePoint;
typedef Point<ExactNumber> ExactPoint;


// Domain declarations
using IntervalDomain = ExactInterval;
using BoxDomain = ExactBox;

// Function declarations
template<class P, class D, class F> class Function;
template<class P, class D=BoxDomain> using ScalarFunction = Function<P,D,IntervalDomain>;
template<class P, class D=BoxDomain> using VectorFunction = Function<P,D,BoxDomain>;
template<class P, class C> using UnivariateFunction = Function<P,IntervalDomain,C>;
template<class P, class C> using MultivariateFunction = Function<P,BoxDomain,C>;

template<class P> using ScalarUnivariateFunction = Function<P,IntervalDomain,IntervalDomain>;
template<class P> using VectorUnivariateFunction = Function<P,IntervalDomain,BoxDomain>;
template<class P> using ScalarMultivariateFunction = Function<P,BoxDomain,IntervalDomain>;
template<class P> using VectorMultivariateFunction = Function<P,BoxDomain,BoxDomain>;

typedef ScalarFunction<ApproximateTag> ApproximateScalarFunction;
typedef ScalarFunction<ValidatedTag> ValidatedScalarFunction;
typedef ScalarFunction<EffectiveTag> EffectiveScalarFunction;
typedef EffectiveScalarFunction RealScalarFunction;

typedef VectorFunction<ApproximateTag> ApproximateVectorFunction;
typedef VectorFunction<ValidatedTag> ValidatedVectorFunction;
typedef VectorFunction<EffectiveTag> EffectiveVectorFunction;
typedef EffectiveVectorFunction RealVectorFunction;

// Function interface declarations
template<class P, class D, class C> class FunctionInterface;
template<class P, class D=BoxDomain> using ScalarFunctionInterface = FunctionInterface<P,D,IntervalDomain>;
template<class P, class D=BoxDomain> using VectorFunctionInterface = FunctionInterface<P,D,BoxDomain>;

typedef ScalarFunctionInterface<ApproximateTag> ApproximateScalarFunctionInterface;
typedef ScalarFunctionInterface<ValidatedTag> ValidatedScalarFunctionInterface;
typedef ScalarFunctionInterface<EffectiveTag> EffectiveScalarFunctionInterface;

typedef VectorFunctionInterface<ApproximateTag> ApproximateVectorFunctionInterface;
typedef VectorFunctionInterface<ValidatedTag> ValidatedVectorFunctionInterface;
typedef VectorFunctionInterface<EffectiveTag> EffectiveVectorFunctionInterface;

// Function model declarations
template<class P> class ScalarFunctionModel;
template<class P> class VectorFunctionModel;

typedef ScalarFunctionModel<ApproximateTag> ApproximateScalarFunctionModel;
typedef ScalarFunctionModel<ValidatedTag> ValidatedScalarFunctionModel;

typedef VectorFunctionModel<ApproximateTag> ApproximateVectorFunctionModel;
typedef VectorFunctionModel<ValidatedTag> ValidatedVectorFunctionModel;


} // namespace Ariadne

#endif
