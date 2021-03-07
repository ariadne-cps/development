/***************************************************************************
 *            utility/identifier.hpp
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

/*! \file utility/identifier.hpp
 *  \brief Strings used as identifiers.
 */

#ifndef ARIADNE_IDENTIFIER_HPP
#define ARIADNE_IDENTIFIER_HPP

#include "string.hpp"

namespace Ariadne {

//! \brief A class representing the name of an object.
//! \details A proxy for a standard string.
//! \sa Constant, Variable, String
class Identifier : public std::string
{
  public:
    Identifier() : std::string() { }
    Identifier(const char* cstr) : std::string(cstr) { }
    //! \brief Construct an identifier from a standard string.
    Identifier(const std::string& str) : std::string(str) { }
};

} // namespace Ariadne

#endif /* ARIADNE_IDENTIFIER_HPP */
