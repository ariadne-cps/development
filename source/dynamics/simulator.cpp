/***************************************************************************
 *            dynamics/simulator.cpp
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

#include "../function/functional.hpp"

#include "../algebra/algebra.hpp"

#include "../config.hpp"

#include "../utility/array.hpp"
#include "../utility/container.hpp"
#include "../utility/tuple.hpp"
#include "../utility/stlio.hpp"
#include "../symbolic/valuation.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/space.hpp"

#include "../function/function.hpp"
#include "../function/formula.hpp"
#include "../function/taylor_model.hpp"

#include "../solvers/runge_kutta_integrator.hpp"

#include "../dynamics/orbit.hpp"
#include "../dynamics/vector_field.hpp"
#include "../dynamics/simulator.hpp"

namespace Ariadne {

template class Orbit<Point<FloatDPApproximation>>;

Simulator::Simulator() : _step_size(0.125_x,dp)
{ }

Void Simulator::set_step_size(double h) {
    this->_step_size=h;
}

inline FloatDPApproximation evaluate(const EffectiveScalarMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }
inline Vector<FloatDPApproximation> evaluate(const EffectiveVectorMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }

auto Simulator::orbit(const VectorField& system, const RealPoint& init_pt, const TerminationType& termination) const
    -> Orbit<ApproximatePointType>
{
    return orbit(system,ApproximatePointType(init_pt,dp),termination);
}

auto Simulator::orbit(const VectorField& system, const ApproximatePointType& init_pt, const TerminationType& termination) const
    -> Orbit<ApproximatePointType>
{
    ARIADNE_LOG_SCOPE_CREATE;

    VectorField::TimeType t(0.0_exact);
    Dyadic h(ExactDouble(this->_step_size.get_d()));
    VectorField::TimeType tmax(termination);

    Orbit<ApproximatePointType> orbit(init_pt);

    EffectiveVectorMultivariateFunction dynamic=system.function();

    RungeKutta4Integrator integrator(this->_step_size.get_d());

    ApproximatePointType point=init_pt;

    while(possibly(t<tmax)) {
        Int old_precision = std::clog.precision();
        ARIADNE_LOG_PRINTLN_AT(1,
                "t=" << std::setw(4) << std::left << t.compute(Effort(0u))
                << " p=" << point
                << std::setprecision(old_precision));

        point = ApproximatePointType(integrator.step(dynamic,point,this->_step_size));
        t += h;
        orbit.insert(t.compute_get(Effort(0),DoublePrecision()).value(),point);
    }

    return orbit;
}



}  // namespace Ariadne

