/***************************************************************************
 *            test_verification_manager.cpp
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

#include "verification/verification_manager.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestVerificationManager {
  public:

    void test_instantiation() {
        VerificationManager::instance();
    }

    void test() {
        ARIADNE_TEST_CALL(test_instantiation());
    }
};

int main() {
    TestVerificationManager().test();
    return ARIADNE_TEST_FAILURES;
}
