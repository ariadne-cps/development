/***************************************************************************
 *            test_linear_programming.cpp
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

#include <iostream>
#include <fstream>

#include "config.hpp"
#include "../test.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "solvers/linear_programming.hpp"


using namespace std;
using namespace Ariadne;

double infd=inf.get_d();

class TestInteriorPointSolver
{
    InteriorPointSolver* optimiser;
  public:
    TestInteriorPointSolver(InteriorPointSolver& o) : optimiser(&o) { }

    Void test() {
        ARIADNE_TEST_CALL(test_validate_feasibility());
        //ARIADNE_TEST_CALL(test_optimization());
        //ARIADNE_TEST_CALL(test_feasibility());
        ARIADNE_TEST_CALL(test_constrained_feasibility());
    }

    Void test_validate_feasibility() {

        // A feasible instance
        Vector<FloatDPValue> xl({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDPValue> xu({inf,inf,inf},dp);
        Matrix<FloatDPValue> A({{1.0_x,0.0_x,1.0_x}, {0.0_x,1.0_x,2.0_x}},dp);
        Vector<FloatDPValue> b({1.0_x,1.0_x},dp);
        Vector<FloatDPValue> x({0.81_pr,0.51_pr,0.21_pr},dp);
        Vector<FloatDPValue> y({-1.0_x,-0.5_x},dp);
        ARIADNE_ASSERT(definitely(optimiser->validate_feasibility(xl,xu,A,b,x,y)));

        b=Vector<FloatDPValue>({1.0_x, -1.0_x},dp);
        y=Vector<FloatDPValue>({-0.5_x, -1.0_x},dp);
        ARIADNE_ASSERT(definitely(!optimiser->validate_feasibility(xl,xu,A,b,x,y)));
    }


    Void test_feasibility() {
        // A feasible instance
        Vector<FloatDPValue> xl({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDPValue> xu({inf,inf,inf},dp);
        Matrix<FloatDPValue> A({{1.0_x,0.0_x,1.0_x},{0.0_x,1.0_x,2.0_x}},dp);
        Vector<FloatDPValue> b({1.0_x,1.0_x},dp);
        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=Vector<FloatDPValue>({1.0_x,-1.0_x},dp);
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));
    }


    Void test_constrained_feasibility() {
        Matrix<FloatDPValue> A({{1.0_x,0.0_x,1.0_x},{0.0_x,1.0_x,2.0_x}},dp);
        Vector<FloatDPValue> b({1.0_x,1.0_x},dp);
        Vector<FloatDPValue> xl({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDPValue> xu({4.0_x,2.0_x,3.0_x},dp);

        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=Vector<FloatDPValue>({1.0_x,-1.0_x},dp);
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));

        xu=Vector<FloatDPValue>({+inf,+inf,3.0_x},dp);
    }


    Void test_optimization() {
        Matrix<FloatDPValue> A({{1.0_x,0.0_x,1.0_x},{0.0_x,1.0_x,2.0_x}},dp);
        Vector<FloatDPValue> b({1.0_x,1.0_x},dp);
        Vector<FloatDPValue> c({1.0_x,0.5_x,-0.75_x},dp);
        Vector<FloatDPValue> xl({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDPValue> xu({+inf,+inf,3.0_x},dp);

        ARIADNE_TEST_PRINT(optimiser->minimise(c,xl,xu,A,b));
    }

};


Int main(Int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));

    InteriorPointSolver interior_point_optimiser;
    TestInteriorPointSolver(interior_point_optimiser).test();

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

