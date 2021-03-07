/***************************************************************************
 *            test_concurrency_manager.cpp
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

#include <thread>
#include "concurrency/concurrency_manager.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestConcurrencyManager {
  public:

    void test_check_nonzero_maximum_concurrency() {
        ARIADNE_TEST_ASSERT(ConcurrencyManager::instance().maximum_concurrency()>0);
    }

    void test_set_concurrency() {
        auto max_concurrency = ConcurrencyManager::instance().maximum_concurrency();
        ConcurrencyManager::instance().set_concurrency(max_concurrency);
        ARIADNE_TEST_EQUALS(ConcurrencyManager::instance().concurrency(),max_concurrency);
        ARIADNE_TEST_FAIL(ConcurrencyManager::instance().set_concurrency(0));
        ARIADNE_TEST_FAIL(ConcurrencyManager::instance().set_concurrency(1+max_concurrency));
    }

    void test() {
        ARIADNE_TEST_CALL(test_check_nonzero_maximum_concurrency());
        ARIADNE_TEST_CALL(test_set_concurrency());
    }
};

int main() {
    TestConcurrencyManager().test();
    return ARIADNE_TEST_FAILURES;
}
