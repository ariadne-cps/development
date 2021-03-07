/***************************************************************************
 *            arch_suite.cpp
 *
 *  Copyright  2020  Luca Geretti
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

#include "PRDE20.hpp"
#include "CVDP20.hpp"
#include "LALO20.hpp"
#include "QUAD20.hpp"
#include "LOVO20.hpp"
#include "SPRE20.hpp"

using namespace Ariadne;

Int main(Int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));
    PRDE20();
    CVDP20();
    LALO20();
    QUAD20();
    LOVO20();
    SPRE20();
}
