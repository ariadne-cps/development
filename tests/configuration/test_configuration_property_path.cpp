/***************************************************************************
 *            test_configuration_property_path.cpp
 *
 *  Copyright  2008-20  Luca Geretti
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

#include "utility/identifier.hpp"
#include "configuration/configuration_property_path.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestConfigurationPropertyPath {
  public:

    void test_construction() {
        ConfigurationPropertyPath p;
        ARIADNE_TEST_EQUALS(p.repr(),"./");
        ConfigurationPropertyPath p2("child");
        ARIADNE_TEST_EQUALS(p2.repr(),"./child/");
        ConfigurationPropertyPath p3(p);
        ARIADNE_TEST_EQUALS(p3.repr(),"./");
    }

    void test_append() {
        ConfigurationPropertyPath p;
        p.append("child1");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/");
        p.append("child2");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/child2/");
    }

    void test_prepend() {
        ConfigurationPropertyPath p;
        p.prepend("child2");
        ARIADNE_TEST_EQUALS(p.repr(),"./child2/");
        p.prepend("child1");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/child2/");
    }

    void test_first_last_subpath() {
        ConfigurationPropertyPath p;
        ARIADNE_TEST_FAIL(p.first());
        ARIADNE_TEST_FAIL(p.last())
        p.append("child1");
        p.append("child2");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/child2/");
        auto sp = p.subpath();
        ARIADNE_TEST_EQUALS(sp.repr(),"./child2/");
        auto f = p.first();
        ARIADNE_TEST_EQUALS(f,"child1");
        auto l = p.last();
        ARIADNE_TEST_EQUALS(l,"child2");
    }

    void test_copy() {
        ConfigurationPropertyPath p1;
        p1.append("child1");
        auto p2 = p1;
        p2.append("child2");
        ARIADNE_TEST_EQUALS(p1.repr(),"./child1/");
        ARIADNE_TEST_EQUALS(p2.repr(),"./child1/child2/");
    }

    void test_less_equal() {
        ConfigurationPropertyPath p1;
        p1.append("child1");
        ConfigurationPropertyPath p2;
        p2.append("child1");
        ARIADNE_TEST_EQUAL(p1,p2);
        p2.append("child2");
        ARIADNE_TEST_ASSERT(p1<p2);
        p2.prepend("child0");
        ARIADNE_TEST_ASSERT(p2<p1);
    }

    void test() {
        ARIADNE_TEST_CALL(test_construction());
        ARIADNE_TEST_CALL(test_append());
        ARIADNE_TEST_CALL(test_prepend());
        ARIADNE_TEST_CALL(test_first_last_subpath());
        ARIADNE_TEST_CALL(test_copy());
        ARIADNE_TEST_CALL(test_less_equal());
    }
};

int main() {
    TestConfigurationPropertyPath().test();
    return ARIADNE_TEST_FAILURES;
}
