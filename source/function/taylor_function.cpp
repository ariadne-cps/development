/***************************************************************************
 *            function/taylor_function.cpp
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

#include "function/functional.hpp"
#include "config.hpp"

#include <iostream>
#include <iomanip>

#include "utility/macros.hpp"
#include "utility/exceptions.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "algebra/multi_index.hpp"
#include "function/polynomial.hpp"
#include "algebra/differential.hpp"
#include "algebra/evaluate.hpp"
#include "function/taylor_model.hpp"

#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/scaled_function_patch.hpp"

#include "function/taylor_function.hpp"

#include "taylor_model.tpl.hpp"
#include "scaled_function_patch.tpl.hpp"
#include "function_mixin.tpl.hpp"

#define VOLATILE ;

namespace Ariadne {

template class ScaledFunctionPatchFactory<ValidatedTaylorModelDP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelDP>,RealVector>;

template class ScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,RealScalar(RealVector)>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,RealScalar(RealVector)>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,RealVector(RealVector)>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,RealVector(RealVector)>;

template class ScaledFunctionPatchFactory<ValidatedTaylorModelMP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelMP>,RealVector>;

template class ScaledFunctionPatch<ValidatedTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelMP>;


FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper<FloatDP>());
}

FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory(Sweeper<FloatDP> const& sweeper) {
    return new TaylorFunctionFactory(sweeper);
}

} // namespace Ariadne
