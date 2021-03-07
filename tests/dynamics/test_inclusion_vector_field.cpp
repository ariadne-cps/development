/***************************************************************************
 *            test_inclusion_vector_field.cpp
 *
 *  Copyright  2008-20  Luca Geretti, Pieter Collins, Sanja Zivanovic
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


#include "geometry/box.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "symbolic/expression_set.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/symbolic_function.hpp"
#include "geometry/function_set.hpp"

#include "dynamics/inclusion_vector_field.hpp"


#include "../test.hpp"

using namespace Ariadne;

class TestInclusionVectorField {
  public:

    Void test_not_formula_function() {
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(3,1);

        EffectiveVectorMultivariateFunction dynamics({x,y});

        BoxDomainType inputs(1,{1,2});

        ARIADNE_TEST_THROWS(InclusionVectorField(dynamics,inputs),NotFormulaFunctionException);
    }

    Void test_missing_input() {
        RealVariable x("x"), y("y"), z("z"), u("u");
        DottedRealAssignments dynamics={dot(x)=y+u,dot(y)=-x+z};
        RealVariableIntervals inputs={-4/100_q<=u<=4/100_q};

        ARIADNE_TEST_THROWS(InclusionVectorField(dynamics,inputs),MissingInputException);
    }

    Void test_unused_input() {
        RealVariable x("x"), y("y"), u("u");
        DottedRealAssignments dynamics={dot(x)=y,dot(y)=-x};
        RealVariableIntervals inputs={-4/100_q<=u<=4/100_q};

        ARIADNE_TEST_THROWS(InclusionVectorField(dynamics,inputs),UnusedInputException);
    }

    Void test_mismatching_coordinates() {
        RealVariable x("x"), y("y");
        DottedRealAssignments dynamics={dot(x)=-y,dot(y)=x};

        EffectiveVectorMultivariateFunction function = make_function(left_hand_sides(dynamics),Vector<RealExpression>(right_hand_sides(dynamics)));

        BoxDomainType inputs(1,{1,2});

        ARIADNE_TEST_THROWS(InclusionVectorField(function,inputs),FunctionArgumentsMismatchException);
    }

    Void test_dynamics_transformation() {
        RealVariable x("x"), y("y"), u("u");
        DottedRealAssignments dynamics={dot(x)=-y,dot(y)=x+2*u};
        RealVariableIntervals inputs={1<=u<=2};

        InclusionVectorField ivf(dynamics,inputs);

        ARIADNE_TEST_PRINT(ivf);

        auto nic = ivf.noise_independent_component();
        const EffectiveVectorFormulaFunction& noise_independent_component = dynamic_cast<const EffectiveVectorFormulaFunction&>(nic.reference());
        ARIADNE_TEST_PRINT(noise_independent_component);

        EffectiveFormula zero = EffectiveFormula::zero();
        EffectiveFormula one = EffectiveFormula::constant(1_dec);
        EffectiveFormula twice = EffectiveFormula::constant(2_dec);
        EffectiveFormula one_point_five = EffectiveFormula::constant(1.5_dec);
        EffectiveFormula xf = EffectiveFormula::coordinate(0);
        EffectiveFormula yf = EffectiveFormula::coordinate(1);

        ARIADNE_TEST_ASSERT(identical(noise_independent_component._formulae,{-yf,xf+twice*one_point_five}));

        auto id = ivf.input_derivatives();
        const EffectiveVectorFormulaFunction& input_derivatives = dynamic_cast<const EffectiveVectorFormulaFunction&>(id[0].reference());
        ARIADNE_TEST_PRINT(input_derivatives);

        ARIADNE_TEST_ASSERT(identical(input_derivatives._formulae,{zero,one}));

        ARIADNE_TEST_EQUAL(ivf.dimension(),2);
        ARIADNE_TEST_EQUAL(ivf.number_of_inputs(),1);
        ARIADNE_TEST_ASSERT(ivf.is_input_affine());
        ARIADNE_TEST_ASSERT(ivf.is_input_additive());
    }

    Void test_nonaffine_field() {
        RealVariable x("x"), y("y"), u1("u1"), u2("u2");
        DottedRealAssignments dynamics={dot(x)=u1*u2*x*(1-y),dot(y)=u2*y*(x-1)};
        RealVariableIntervals inputs={2.99_dec<=u1<=3.01_dec,0.99_dec<=u2<=1.01_dec};

        InclusionVectorField ivf(dynamics,inputs);

        ARIADNE_TEST_PRINT(ivf);
        ARIADNE_TEST_PRINT(ivf.noise_independent_component());
        ARIADNE_TEST_PRINT(ivf.input_derivatives());

        ARIADNE_TEST_EQUAL(ivf.dimension(),2);
        ARIADNE_TEST_EQUAL(ivf.number_of_inputs(),2);
        ARIADNE_TEST_ASSERT(not ivf.is_input_affine());
        ARIADNE_TEST_ASSERT(not ivf.is_input_additive());
    }

    Void test_nonadditive_field() {
        RealVariable x("x"), y("y"), u1("u1"), u2("u2");
        DottedRealAssignments dynamics={dot(x)=u1*x*(1-y),dot(y)=u2*y*(x-1)};
        RealVariableIntervals inputs={2.99_dec<=u1<=3.01_dec,0.99_dec<=u2<=1.01_dec};

        InclusionVectorField ivf(dynamics,inputs);

        ARIADNE_TEST_PRINT(ivf);
        ARIADNE_TEST_PRINT(ivf.noise_independent_component());
        ARIADNE_TEST_PRINT(ivf.input_derivatives());

        ARIADNE_TEST_EQUAL(ivf.dimension(),2);
        ARIADNE_TEST_EQUAL(ivf.number_of_inputs(),2);
        ARIADNE_TEST_ASSERT(ivf.is_input_affine());
        ARIADNE_TEST_ASSERT(not ivf.is_input_additive());
    }

    Void test_additive_field() {
        RealVariable x("x"), y("y"), z("z"), u("u");
        DottedRealAssignments dynamics={dot(x)=-y-z,dot(y)=x+y*0.1_dec,dot(z)=z*(x-6)+u};
        RealVariableIntervals inputs={0.099_dec<=u<=0.101_dec};

        InclusionVectorField ivf(dynamics,inputs);

        ARIADNE_TEST_PRINT(ivf);
        ARIADNE_TEST_PRINT(ivf.noise_independent_component());
        ARIADNE_TEST_PRINT(ivf.input_derivatives());

        ARIADNE_TEST_EQUAL(ivf.dimension(),3);
        ARIADNE_TEST_EQUAL(ivf.number_of_inputs(),1);
        ARIADNE_TEST_ASSERT(ivf.is_input_affine());
        ARIADNE_TEST_ASSERT(ivf.is_input_additive());
    }

    Void test() {
        ARIADNE_TEST_CALL(test_not_formula_function());
        ARIADNE_TEST_CALL(test_missing_input());
        ARIADNE_TEST_CALL(test_unused_input());
        ARIADNE_TEST_CALL(test_mismatching_coordinates());
        ARIADNE_TEST_CALL(test_dynamics_transformation());
        ARIADNE_TEST_CALL(test_nonaffine_field());
        ARIADNE_TEST_CALL(test_nonadditive_field());
        ARIADNE_TEST_CALL(test_additive_field());
    }
};

int main() {

    TestInclusionVectorField().test();
    return ARIADNE_TEST_FAILURES;
}
