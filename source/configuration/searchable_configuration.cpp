/***************************************************************************
 *            concurrency/searchable_configuration.cpp
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

#include "searchable_configuration.hpp"
#include "configuration_search_space.hpp"

namespace Ariadne {

SearchableConfiguration::SearchableConfiguration(SearchableConfiguration const& c) {
    for (auto p : c.properties()) _properties.insert(Pair<Identifier,SharedPointer<ConfigurationPropertyInterface>>(p.first,SharedPointer<ConfigurationPropertyInterface>(p.second->clone())));
}

SearchableConfiguration& SearchableConfiguration::operator=(SearchableConfiguration const& c) {
    _properties.clear();
    for (auto p : c.properties()) _properties.insert(Pair<Identifier,SharedPointer<ConfigurationPropertyInterface>>(p.first,SharedPointer<ConfigurationPropertyInterface>(p.second->clone())));
    return *this;
}

Map<Identifier,SharedPointer<ConfigurationPropertyInterface>>& SearchableConfiguration::properties() {
    return _properties;
}

Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> const& SearchableConfiguration::properties() const {
    return _properties;
}

void SearchableConfiguration::add_property(Identifier const& name, ConfigurationPropertyInterface const& property) {
    _properties.insert(Pair<Identifier,SharedPointer<ConfigurationPropertyInterface>>({name,SharedPointer<ConfigurationPropertyInterface>(property.clone())}));
}

OutputStream& SearchableConfiguration::_write(OutputStream& os) const {
    os << "(\n";
    auto iter = _properties.begin(); SizeType i=0;
    while (i<_properties.size()-1) {
        os << iter->first << " = " << *iter->second << ",\n";
        ++iter; ++i;
    }
    os << iter->first << " = " << *iter->second << ")"; return os;
}

Bool SearchableConfiguration::is_singleton() const {
    for (auto p : _properties) {
        auto int_values = p.second->integer_values();
        for (auto p_int : int_values) if (p_int.second.size() > 1) return false;
    }
    return true;
}

ConfigurationSearchSpace SearchableConfiguration::search_space() const {
    Set<ConfigurationSearchParameter> result;
    for (auto p : _properties) {
        auto integer_values = p.second->integer_values();
        for (auto p_int : integer_values) {
            if (p_int.second.size() > 1) {
                ConfigurationPropertyPath path(p_int.first);
                path.prepend(p.first);
                result.insert(ConfigurationSearchParameter(path, p.second->is_metric(p_int.first), p_int.second));
            }
        }
    }
    ARIADNE_ASSERT_MSG(not result.empty(),"The search space is empty.");
    return result;
}

} // namespace Ariadne
