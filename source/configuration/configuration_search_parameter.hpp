/***************************************************************************
 *            configuration/configuration_search_parameter.hpp
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

/*! \file configuration/configuration_search_parameter.hpp
 *  \brief Classes for handling configuration search parameters.
 */

#ifndef ARIADNE_CONFIGURATION_SEARCH_PARAMETER_HPP
#define ARIADNE_CONFIGURATION_SEARCH_PARAMETER_HPP

#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/macros.hpp"
#include "utility/identifier.hpp"
#include "configuration_property_path.hpp"

namespace Ariadne {

class ConfigurationSearchParameter {
  public:
    ConfigurationSearchParameter(ConfigurationPropertyPath const& path, Bool is_metric, List<int> const& values);
    ConfigurationPropertyPath const& path() const;
    //! \brief Admissible values
    List<int> const& values() const;
    //! \brief Whether the parameter should shift to adjacent values instead of hopping between values
    Bool is_metric() const;
    //! \brief Generate a random value, useful for the initial value
    int random_value() const;
    //! \brief Randomly get the result from shifting the given \a value
    int shifted_value_from(int value) const;

    Bool operator==(ConfigurationSearchParameter const& p) const;
    Bool operator<(ConfigurationSearchParameter const& p) const;

    friend OutputStream& operator<<(OutputStream& os, ConfigurationSearchParameter const& parameter);

  private:
    const ConfigurationPropertyPath _path;
    const Bool _is_metric;
    const List<int> _values;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_SEARCH_PARAMETER_HPP
