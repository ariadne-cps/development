/***************************************************************************
 *            concurrency/configuration_property.cpp
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

#include <ostream>
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "configuration_property.hpp"
#include "configuration_property.tpl.hpp"

namespace Ariadne {

BooleanConfigurationProperty::BooleanConfigurationProperty()
    : ConfigurationPropertyBase(false), _is_single(false)
{ }

BooleanConfigurationProperty::BooleanConfigurationProperty(Bool const& value)
    : ConfigurationPropertyBase(true), _is_single(true), _value(value)
{ }

ConfigurationPropertyInterface* BooleanConfigurationProperty::at(ConfigurationPropertyPath const& path) {
    ARIADNE_ASSERT_MSG(path.is_root(),"The path " << path << " is not a root but a boolean property can't have configurable objects below.");
    return this;
}

Bool const& BooleanConfigurationProperty::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return _value;
}

Bool BooleanConfigurationProperty::is_single() const {
    return _is_single;
}

Bool BooleanConfigurationProperty::is_metric(ConfigurationPropertyPath const& path) const {
    ARIADNE_PRECONDITION(path.is_root());
    return false;
}

Bool BooleanConfigurationProperty::is_configurable() const {
    return false;
}

SizeType BooleanConfigurationProperty::cardinality() const {
    if (_is_single) return 1;
    else if (this->is_specified()) return 2;
    else return 0;
}

List<int> BooleanConfigurationProperty::local_integer_values() const {
    List<int> result;
    if (_is_single) result.push_back(_value);
    else if (is_specified()) { result.push_back(0); result.push_back(1); }
    return result;
}

void BooleanConfigurationProperty::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    ARIADNE_PRECONDITION(path.is_root());
    local_set_single(integer_value);
}

void BooleanConfigurationProperty::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not _is_single);
    ARIADNE_PRECONDITION(integer_value == 0 or integer_value == 1);
    _is_single = true;
    if (integer_value == 1) _value = true;
    else _value = false;
}

ConfigurationPropertyInterface* BooleanConfigurationProperty::clone() const {
    return new BooleanConfigurationProperty(*this);
}

void BooleanConfigurationProperty::set_both() {
    set_specified();
    _is_single = false;
}

void BooleanConfigurationProperty::set(Bool const& value) {
    set_specified();
    _is_single = true;
    _value=value;
}

List<SharedPointer<Bool>> BooleanConfigurationProperty::values() const {
    List<SharedPointer<Bool>> result;
    if (is_specified()) {
        if (_is_single) result.append(SharedPointer<Bool>(new Bool(_value)));
        else {
            result.append(SharedPointer<Bool>(new Bool(true)));
            result.append(SharedPointer<Bool>(new Bool(false)));
        }
    }
    return result;
}



} // namespace Ariadne
