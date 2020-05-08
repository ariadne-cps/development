/***************************************************************************
 *            solvers/bounder.cpp
 *
 *  Copyright  2018-20  Luca Geretti
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

#include "bounder.hpp"
#include "../function/formula.hpp"
#include "../function/taylor_model.hpp"

namespace Ariadne {

Pair<StepSizeType,UpperBoxType> BounderBase::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const {
    return this->compute(f,D,BoxDomainType(0u),hsug);
}

Pair<StepSizeType,UpperBoxType> BounderBase::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const {
    return this->compute(f,D,t,BoxDomainType(0u),hsug);
}


Pair<StepSizeType,UpperBoxType> EulerBounder::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const {
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+A.dimension());
    return this->_compute(f,D,0,A,hsug);
}

Pair<StepSizeType,UpperBoxType> EulerBounder::compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const {
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+1u+A.dimension());
    return this->_compute(f,D,t,A,hsug);
}

Pair<StepSizeType,UpperBoxType> EulerBounder::_compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const {
    const PositiveFloatDPValue BOX_RADIUS_WIDENING=cast_positive(0.25_exact);
    const PositiveFloatDPValue NO_WIDENING=cast_positive(1.0_exact);
    const PositiveFloatDPValue INITIAL_STARTING_WIDENING=cast_positive(2.0_exact);
    const PositiveFloatDPValue INITIAL_REFINING_WIDENING=cast_positive(1.125_exact);
    const PositiveFloatDPValue LIPSCHITZ_TOLERANCE=cast_positive(1.0_exact);
    const StepSizeType MINIMUM_STEP_SIZE(1,20u);
    const Nat EXPANSION_STEPS=4;
    const Nat H_REFINEMENT_STEPS=30;
    const Nat REFINEMENT_STEPS=4;

    StepSizeType h=hsug;

    FloatDPUpperBound lipschitz = norm(f.jacobian(Vector<FloatDPBounds>(cast_singleton(product(D,to_time_bounds(t,t+h),A))))).upper();
    StepSizeType hlip = static_cast<StepSizeType>(cast_exact(LIPSCHITZ_TOLERANCE/lipschitz));
    h=min(hlip,h);

    IntervalDomainType T = to_time_bounds(t,t+h);

    UpperBoxType B;
    Nat h_refinement_step = 0;
    Dyadic result_h(0);
    UpperBoxType result_B = B;
    //std::cout << "Starting from " << h << std::endl;
    while(true) {
        //std::cout << "Check " << h << std::endl;
        B=this->_formula(f,D,T,A,D,BOX_RADIUS_WIDENING,INITIAL_STARTING_WIDENING);
        Bool success=false;
        for(Nat i=0; i<EXPANSION_STEPS; ++i) {
            UpperBoxType Br=this->_refinement(f,D,T,A,B);
            if(not definitely(is_bounded(Br))) {
                success=false;
                break;
            } else if(refines(Br,B)) {
                B=Br;
                success=true;
                break;
            } else {
                B=this->_formula(f,D,T,A,B, NO_WIDENING,INITIAL_REFINING_WIDENING);
            }
        }
        if(!success) {
            //std::cout << "Failure" << std::endl;
            h=hlf(h);
            if (h < MINIMUM_STEP_SIZE)
                ARIADNE_THROW(BoundingNotFoundException,"EulerBounder::_compute","The step size is lower than the minimum (" << MINIMUM_STEP_SIZE << ") allowed, bounding could not be found.");

            T = to_time_bounds(t,t+h);
        } else {
            //std::cout << "Success" << std::endl;
            if (h > result_h) {
                result_h = h;
                result_B = B;
            }
            ++h_refinement_step;
            if (h_refinement_step >= H_REFINEMENT_STEPS) break;
            else {
                h = hlf(3*h);
                T = to_time_bounds(t,t+h);
            }
        }
    }

    for(Nat i=0; i<REFINEMENT_STEPS; ++i) {
        result_B = this->_refinement(f,D,T,A,result_B);
    }

    //Dyadic::set_default_writer(new DecimalWriter());
    //std::cout << "using " << result_h << ", " << result_B << std::endl;
    return std::make_pair(result_h,result_B);
}

UpperBoxType EulerBounder::_refinement(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const {
    const PositiveFloatDPValue NO_WIDENING=cast_positive(1.0_exact);
    return _formula(f,D,T,A,B, NO_WIDENING,NO_WIDENING);
}

UpperBoxType EulerBounder::_formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B, PositiveFloatDPValue INITIAL_BOX_WIDENING, PositiveFloatDPValue VECTOR_WIDENING) const {
    UpperIntervalType const& rT=reinterpret_cast<UpperIntervalType const&>(T);
    UpperBoxType const& rA=reinterpret_cast<UpperBoxType const&>(A);
    UpperIntervalType const rH=rT-T.lower();

    const bool is_autonomous = (f.argument_size() == D.dimension()+A.dimension());
    UpperBoxType dom = is_autonomous ? product(B,rA) : product(B,rT,rA);

    UpperBoxType wD = (INITIAL_BOX_WIDENING!=1) ? D : D + INITIAL_BOX_WIDENING*(D-D.midpoint());

    return wD+(VECTOR_WIDENING*rH)*cast_vector(apply(f,dom));
}

} // namespace Ariadne;
