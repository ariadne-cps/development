/***************************************************************************
 *            numeric/concepts.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include "concepts.hpp"
#include "archetypes.hpp"

#include "builtin.hpp"
#include "dyadic.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float_bounds.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "float_approximation.hpp"

namespace Ariadne {

static_assert(IsRoundedLatticeField<FloatDP>);
static_assert(IsRoundedLatticeField<FloatMP>);
static_assert(IsRoundedLatticeField<RoundedArchetype>);



template<> Nat Approximation<RoundedArchetype>::output_places=0;

//template class Bounds<RoundedArchetype>;
template class LowerBound<RoundedArchetype>;
template class UpperBound<RoundedArchetype>;
template class Approximation<RoundedArchetype>;

} // namespace Ariadne
