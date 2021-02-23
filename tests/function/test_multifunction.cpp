/***************************************************************************
 *            test_multifunction.cpp
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.hpp"
#include "algebra/sweeper.hpp"
#include "function/multifunction.hpp"
#include "function/taylor_multifunction.hpp"
#include "numeric/numeric.hpp"

#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

class TestMultifunction
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_evaluate();
    Void test_taylor_evaluate();
};

Void TestMultifunction::test()
{
    Dyadic::set_default_writer(DecimalWriter());
    ARIADNE_TEST_CALL(test_evaluate());
    ARIADNE_TEST_CALL(test_taylor_evaluate());
}

Void TestMultifunction::test_concept()
{
}

Void TestMultifunction::test_evaluate()
{
    RealVariable x0("x0"), x1("x1"), p("p");

    BoxDomainType dom({{-1,3},{-1,3}});
    BoxDomainType prms({{-1,1}});;

    EffectiveVectorMultivariateFunction f=make_function({x0,x1,p},{1.5_dy-x0*x0-0.375_dy*x1+p,x0});

    Sweeper<FloatDP> swp;

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionPatch,mvfp,(f,prms));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunction,mvf,(mvfp));

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionModel<DP>,mvfmdp,(dom,f,prms,swp));;
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionModel<DP>,mvfdp,(mvfmdp));;

    Vector<ValidatedNumber> v(2,1);
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedLocatedSet<RealVector>,res,mvfp(v));
    ARIADNE_TEST_ASSIGN(res,mvf(v));
    ARIADNE_TEST_ASSIGN(res,mvfmdp(v));
    ARIADNE_TEST_ASSIGN(res,mvfdp(v));
}

Void TestMultifunction::test_taylor_evaluate()
{
    using ScalarIntervalFunctionModel = ValidatedIntervalTaylorFunctionModel<FloatMP>;
    using VectorIntervalFunctionModel = ValidatedVectorIntervalTaylorFunctionModel<FloatMP>;
    using Ivl = Interval<FloatMPUpperBound>;

    MP pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1}});
//    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    ScalarIntervalFunctionModel x0=ScalarIntervalFunctionModel::coordinate(dom,0,swp);
    ScalarIntervalFunctionModel x1=ScalarIntervalFunctionModel::coordinate(dom,1,swp);

//    Vector<ValidatedNumber> u({0.75_x,2.5_x});
    Vector<ValidatedNumber> u({0.75_x,-0.5_x});
    Vector<FloatMPBounds> v(u,pr);
    Vector<FloatDPBounds> w(u,dp);

    ScalarIntervalFunctionModel sf=Ivl(3,5,pr)*x0*x1+Ivl(7,11,pr);
    std::cout << "x0=" << x0 << "\n";
    std::cout << "sf=" << sf << "\n";
    Interval<FloatMPUpperBound> fv=sf(v);

    std::cout << "sf(u)=" << sf(u) << "\n";
    std::cout << "sf(v)=" << sf(v) << "\n";
    std::cout << "sf(w)=" << sf(w) << "\n";

    VectorIntervalFunctionModel vf({sf,x0});
    std::cout << "vf=" << vf << "\n";
    Box<Interval<FloatMPUpperBound>> vfv=vf(v);

    std::cout << "vf(u)=" << vf(u) << "\n";
    std::cout << "vf(v)=" << vf(v) << "\n";
    std::cout << "vf(w)=" << vf(w) << "\n";
}


Int main() {
    TestMultifunction().test();

    return ARIADNE_TEST_FAILURES;
}

