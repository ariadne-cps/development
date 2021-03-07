/***************************************************************************
 *            configuration/configuration_search_space.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file configuration/configuration_search_space.hpp
 *  \brief Class for handling a search space of configuration properties.
 */

#ifndef ARIADNE_CONFIGURATION_SEARCH_SPACE_HPP
#define ARIADNE_CONFIGURATION_SEARCH_SPACE_HPP

#include "configuration_search_parameter.hpp"

namespace Ariadne {

class ConfigurationSearchPoint;
class Real;
template<class R> class Variable;

using ParameterBindingsMap = Map<ConfigurationPropertyPath,int>;

class ConfigurationSearchSpace {
  public:
    ConfigurationSearchSpace(Set<ConfigurationSearchParameter> const& parameters);

    ConfigurationSearchPoint make_point(ParameterBindingsMap const& bindings) const;
    ConfigurationSearchPoint initial_point() const;

    List<ConfigurationSearchParameter> const& parameters() const;

    //! \brief The total number of points identified by the space
    SizeType total_points() const;
    //! \brief The number of parameters in the space
    SizeType dimension() const;
    //! \brief The index of the given parameter in the ordered space
    SizeType index(ConfigurationSearchParameter const& p) const;
    //! \brief The index of the given parameter identifier in the ordered space
    SizeType index(ConfigurationPropertyPath const& name) const;
    //! \brief The parameter corresponding to the path \a path
    ConfigurationSearchParameter const& parameter(ConfigurationPropertyPath const& path) const;

    ConfigurationSearchSpace* clone() const;

    friend OutputStream& operator<<(OutputStream& os, ConfigurationSearchSpace const& space);

  private:
    const List<ConfigurationSearchParameter> _parameters;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_SEARCH_SPACE_HPP
