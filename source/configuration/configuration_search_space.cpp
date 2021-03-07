/***************************************************************************
 *            configuration/configuration_search_space.cpp
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

#include "configuration/configuration_property_path.hpp"
#include "configuration_search_point.hpp"
#include "configuration_search_space.hpp"

namespace Ariadne {

ConfigurationSearchSpace::ConfigurationSearchSpace(Set<ConfigurationSearchParameter> const& parameters)
        : _parameters(parameters) {
    ARIADNE_PRECONDITION(not _parameters.empty());
}

ConfigurationSearchPoint ConfigurationSearchSpace::make_point(ParameterBindingsMap const& bindings) const {
    ARIADNE_PRECONDITION(bindings.size() == this->dimension())
    ParameterBindingsMap pb;
    for (auto p : _parameters) {
        int v = bindings.find(p.path())->second;
        pb.insert(Pair<ConfigurationPropertyPath,int>(p.path(), v));
    }
    return ConfigurationSearchPoint(*this, pb);
}

ConfigurationSearchPoint ConfigurationSearchSpace::initial_point() const {
    ParameterBindingsMap pb;
    for (auto p : _parameters) {
        pb.insert(Pair<ConfigurationPropertyPath,int>(p.path(), p.random_value()));
    }
    return ConfigurationSearchPoint(*this, pb);
}

SizeType ConfigurationSearchSpace::index(ConfigurationSearchParameter const& p) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i) == p) return i;
    ARIADNE_FAIL_MSG("Task parameter '" << p << "' not found in the space.");
}

SizeType ConfigurationSearchSpace::index(ConfigurationPropertyPath const& path) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i).path() == path) return i;
    ARIADNE_FAIL_MSG("Task parameter with path '" << path << "' not found in the space.");
}

ConfigurationSearchParameter const& ConfigurationSearchSpace::parameter(ConfigurationPropertyPath const& path) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i).path() == path) return _parameters.at(i);
    ARIADNE_FAIL_MSG("Task parameter with path '" << path << "' not found in the space.");
}

List<ConfigurationSearchParameter> const& ConfigurationSearchSpace::parameters() const {
    return _parameters;
}

SizeType ConfigurationSearchSpace::total_points() const {
    SizeType result = 1;
    for (auto p : _parameters) result *= p.values().size();
    return result;
}

SizeType ConfigurationSearchSpace::dimension() const {
    return _parameters.size();
}

ConfigurationSearchSpace* ConfigurationSearchSpace::clone() const {
    return new ConfigurationSearchSpace(*this);
}

OutputStream& operator<<(OutputStream& os, ConfigurationSearchSpace const& space) {
    return os << space._parameters;
}

} // namespace Ariadne
