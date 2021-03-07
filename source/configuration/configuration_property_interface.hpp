/***************************************************************************
 *            concurrency/configuration_property_interface.hpp
 *
 *  Copyright  2011-20  Luca Geretti
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

/*! \file concurrency/configuration_property_interface.hpp
 *  \brief Interface for properties of a configuration.
 */

#ifndef ARIADNE_CONFIGURATION_PROPERTY_INTERFACE_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_INTERFACE_HPP

#include "utility/writable.hpp"
#include "utility/container.hpp"

namespace Ariadne {

class ConfigurationPropertyPath;

class ConfigurationPropertyInterface : public WritableInterface {
  public:
    //! \brief If only one value is specified
    virtual Bool is_single() const = 0;
    //! \brief If values are specified at all
    virtual Bool is_specified() const = 0;
    //! \brief If the property class at the \a path is metric
    virtual Bool is_metric(ConfigurationPropertyPath const& path) const = 0;
    //! \brief If the property object is a configurable itself
    virtual Bool is_configurable() const = 0;
    //! \brief The number of values stored for the property
    //! \details Returns 1 if single, 0 if not specified.
    virtual SizeType cardinality() const = 0;
    //! \brief Set to a single value a given path, starting from this property
    //! \details Supports the storage of objects that are Configurable themselves
    virtual void set_single(ConfigurationPropertyPath const& path, int integer_value) = 0;
    //! \brief The integer values for each property including the current one
    //! \details Supports the storage of objects that are Configurable themselves
    virtual Map<ConfigurationPropertyPath,List<int>> integer_values() const = 0;
    //! \brief Retrieve a pointer to the property at the given \a path
    virtual ConfigurationPropertyInterface* at(ConfigurationPropertyPath const& path) = 0;

    virtual ConfigurationPropertyInterface* clone() const = 0;
    virtual ~ConfigurationPropertyInterface() = default;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_INTERFACE_HPP
