/***************************************************************************
 *            configuration/searchable_configuration.hpp
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

/*! \file configuration/searchable_configuration.hpp
 *  \brief Classes for configuration and related helpers.
 */

#ifndef ARIADNE_SEARCHABLE_CONFIGURATION_HPP
#define ARIADNE_SEARCHABLE_CONFIGURATION_HPP

#include <ostream>
#include <type_traits>
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "utility/identifier.hpp"
#include "configuration_interface.hpp"
#include "configuration_property_interface.hpp"
#include "configuration_property_path.hpp"

namespace Ariadne {

class ConfigurationSearchSpace;

//! \brief Extension of ConfigurationInterface to deal with search in the properties space
class SearchableConfiguration : public ConfigurationInterface {
  public:
    SearchableConfiguration() = default;
    SearchableConfiguration(SearchableConfiguration const& c);
    SearchableConfiguration& operator=(SearchableConfiguration const& c);
    virtual ~SearchableConfiguration() = default;

    //! \brief Construct a search space from the current configuration
    ConfigurationSearchSpace search_space() const;

    //! \brief If the configuration is made of single values
    Bool is_singleton() const;

    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>>& properties();
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> const& properties() const;

    //! \brief Accessors for get and set of a property identified by a path \a path with type \a P
    //! \details Used in practice to get/set properties for verification
    template<class P> P& at(ConfigurationPropertyPath const& path) {
        auto prop_ptr = _properties.find(path.first());
        ARIADNE_ASSERT_MSG(prop_ptr != _properties.end(),"The property '" << path.first() << "' was not found in the configuration.");
        auto p_ptr = dynamic_cast<P*>(prop_ptr->second->at(path.subpath()));
        ARIADNE_ASSERT_MSG(p_ptr != nullptr, "Invalid property cast, check the property class with respect to the configuration created.")
        return *p_ptr;
    }
    template<class P> P const& at(ConfigurationPropertyPath const& path) const {
        auto prop_ptr = _properties.find(path.first());
        ARIADNE_ASSERT_MSG(prop_ptr != _properties.end(),"The property '" << path.first() << "' was not found in the configuration.");
        auto p_ptr = dynamic_cast<P*>(prop_ptr->second->at(path.subpath()));
        ARIADNE_ASSERT_MSG(p_ptr != nullptr, "Invalid property cast, check the property class with respect to the configuration created.")
        return *p_ptr;
    }
    //! \brief Easier accessors for get and set starting from an \a identifier, to be used for Configuration class accessors
    template<class P> P& at(Identifier const& identifier) {
        return at<P>(ConfigurationPropertyPath(identifier));
    }
    template<class P> P const& at(Identifier const& identifier) const {
        return at<P>(ConfigurationPropertyPath(identifier));
    }

    //! \brief Add a property to the configuration
    void add_property(Identifier const& name, ConfigurationPropertyInterface const& property);

    OutputStream& _write(OutputStream& os) const override;
  private:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> _properties;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
