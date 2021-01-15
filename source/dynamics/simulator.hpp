/***************************************************************************
 *            dynamics/simulator.hpp
 *
 *  Copyright  2009-20  Pieter Collins
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

/*! \file dynamics/simulator.hpp
 *  \brief Simulator for continuous systems.
 */

#ifndef ARIADNE_SIMULATOR_HPP
#define ARIADNE_SIMULATOR_HPP

#include "../output/logging.hpp"
#include "../numeric/float.decl.hpp"
#include "../numeric/floatdp.hpp"
#include "../geometry/point.hpp"

namespace Ariadne {

class VectorField;

template<class T> class Orbit;

/*! \brief A class for computing the evolution of a continuous system.
 */
class Simulator
{
  public:
    typedef Point<FloatDPApproximation> ApproximatePointType;
    typedef ApproximatePointType EnclosureType;
    typedef Real TerminationType;
  private:
    FloatDPApproximation _step_size;
  public:

    //! \brief Default constructor.
    Simulator();

    Void set_step_size(double h);

    //!@{
    //! \name Evolution using points.
    //! \brief Compute an approximation to the orbit set.
    Orbit<ApproximatePointType> orbit(VectorField const& system, ApproximatePointType const& initial_point, TerminationType const& termination) const;
    Orbit<ApproximatePointType> orbit(VectorField const& system, RealPoint const& initial_point, TerminationType const& termination) const;
};



} // namespace Ariadne

#endif // ARIADNE_SIMULATOR_HPP
