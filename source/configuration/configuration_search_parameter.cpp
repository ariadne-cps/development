/***************************************************************************
 *            configuration/configuration_search_parameter.cpp
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

#include "configuration_search_parameter.hpp"

namespace Ariadne {

ConfigurationSearchParameter::ConfigurationSearchParameter(ConfigurationPropertyPath const& path, Bool is_metric, List<int> const& values) :
    _path(path), _is_metric(is_metric), _values(values) {
    ARIADNE_PRECONDITION(values.size()>1);
}

ConfigurationPropertyPath const& ConfigurationSearchParameter::path() const {
    return _path;
}

List<int> const& ConfigurationSearchParameter::values() const {
    return _values;
}

Bool ConfigurationSearchParameter::is_metric() const {
    return _is_metric;
}

int ConfigurationSearchParameter::random_value() const {
    return _values[(SizeType)rand() % _values.size()];
}

int ConfigurationSearchParameter::shifted_value_from(int value) const {
    SizeType num_values = _values.size();
    if (_is_metric) {
        if (value == _values[0]) return value+1;
        if (value == _values[num_values-1]) return value-1;
        if ((SizeType)rand() % 2 == 0) return value+1;
        else return value-1;
    } else {
        SizeType rand_value = (SizeType)rand() % (num_values-1);
        if (_values[rand_value] == value) return _values[num_values-1];
        else return _values[rand_value];
    }
}

Bool ConfigurationSearchParameter::operator==(ConfigurationSearchParameter const& p) const {
    return path() == p.path();
}

Bool ConfigurationSearchParameter::operator<(ConfigurationSearchParameter const& p) const {
    return path() < p.path();
}

OutputStream& operator<<(OutputStream& os, ConfigurationSearchParameter const& p) {
    os << "{'" << p._path << "', is_metric=" << p._is_metric << ", values=";
    if (p._is_metric) os << "[" << p._values[0] << ":" << p._values[p._values.size()-1] << "]";
    else os << p._values;
    return os << "}";
}

} // namespace Ariadne
