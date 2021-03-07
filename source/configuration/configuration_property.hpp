/***************************************************************************
 *            concurrency/configuration_property.hpp
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

/*! \file concurrency/configuration_property.hpp
 *  \brief Classes for properties of a configuration.
 */

#ifndef ARIADNE_CONFIGURATION_PROPERTY_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_HPP

#include <ostream>
#include <type_traits>
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "configuration_interface.hpp"
#include "configuration_property_interface.hpp"
#include "configuration_search_space_converter.hpp"

namespace Ariadne {

template<class T> class ConfigurationPropertyBase : public ConfigurationPropertyInterface {
  protected:
    ConfigurationPropertyBase(Bool const& is_specified);
    void set_specified();
    virtual List<SharedPointer<T>> values() const = 0;
    virtual void local_set_single(int integer_value) = 0;
    virtual List<int> local_integer_values() const = 0;
  public:
    virtual T const& get() const = 0;
    virtual void set(T const& value) = 0;

    Bool is_specified() const override;
    Map<ConfigurationPropertyPath,List<int>> integer_values() const override;

    //! \brief Supplies the values from the property, empty if not specified, the lower/upper bounds if a range
    OutputStream& _write(OutputStream& os) const override;
  private:
    Bool _is_specified;
};

//! \brief A property for a boolean value.
class BooleanConfigurationProperty final : public ConfigurationPropertyBase<Bool> {
  public:
    BooleanConfigurationProperty();
    BooleanConfigurationProperty(Bool const& value);

    Bool is_single() const override;
    Bool is_metric(ConfigurationPropertyPath const& path) const override;
    Bool is_configurable() const override;
    SizeType cardinality() const override;

    ConfigurationPropertyInterface* clone() const override;

    ConfigurationPropertyInterface* at(ConfigurationPropertyPath const& path) override;

    Bool const& get() const override;
    void set(Bool const& value) override;
    void set_both(); //! \brief Set to both true and false
    void set_single(ConfigurationPropertyPath const& path, int integer_value) override;
  protected:
    void local_set_single(int integer_value) override;
    List<int> local_integer_values() const override;
    List<SharedPointer<Bool>> values() const override;
  private:
    Bool _is_single;
    Bool _value;
};

//! \brief A range configuration property offers a range of values with a distance metric
//! \details This property needs a converter to decide how to distribute the integer values in the search space.
template<class T> class RangeConfigurationProperty final : public ConfigurationPropertyBase<T> {
  public:
    RangeConfigurationProperty(ConfigurationSearchSpaceConverterInterface<T> const& converter = LinearSearchSpaceConverter<T>());
    RangeConfigurationProperty(T const& lower, T const& upper, ConfigurationSearchSpaceConverterInterface<T> const& converter = LinearSearchSpaceConverter<T>());
    RangeConfigurationProperty(T const& value, ConfigurationSearchSpaceConverterInterface<T> const& converter = LinearSearchSpaceConverter<T>());

    Bool is_single() const override;
    Bool is_metric(ConfigurationPropertyPath const& path) const override;
    Bool is_configurable() const override;
    SizeType cardinality() const override;

    ConfigurationPropertyInterface* clone() const override;

    ConfigurationPropertyInterface* at(ConfigurationPropertyPath const& path) override;

    T const& get() const override;
    void set(T const& lower, T const& upper);
    //! \brief Set a single value
    //! \details An unbounded single value is accepted
    void set(T const& value) override;
    void set_single(ConfigurationPropertyPath const& path, int integer_value) override;
  protected:
    void local_set_single(int integer_value) override;
    List<int> local_integer_values() const override;
    List<SharedPointer<T>> values() const override;
  private:
    T _lower;
    T _upper;
    SharedPointer<ConfigurationSearchSpaceConverterInterface<T>> const _converter;
};

//! \brief A property that specifies distinct values from an enum
template<class T> class EnumConfigurationProperty final : public ConfigurationPropertyBase<T> {
public:
    EnumConfigurationProperty();
    EnumConfigurationProperty(Set<T> const& values);
    EnumConfigurationProperty(T const& value);

    Bool is_single() const override;
        Bool is_metric(ConfigurationPropertyPath const& path) const override;
    Bool is_configurable() const override;
    SizeType cardinality() const override;

    ConfigurationPropertyInterface* clone() const override;

    ConfigurationPropertyInterface* at(ConfigurationPropertyPath const& path) override;

    T const& get() const override;
    void set(T const& value) override;
    void set(Set<T> const& values);
    void set_single(ConfigurationPropertyPath const& path, int integer_value) override;
protected:
    void local_set_single(int integer_value) override;
    List<int> local_integer_values() const override;
    List<SharedPointer<T>> values() const override;
private:
    Set<T> _values;
};

//! \brief A property that specifies a set of distinct values from handle class \a T
//! \details This can be used either for an enum or for distinct objects of a class or handle class
template<class T> class HandleListConfigurationProperty final : public ConfigurationPropertyBase<T> {
public:
    HandleListConfigurationProperty();
    HandleListConfigurationProperty(List<T> const& values);
    HandleListConfigurationProperty(T const& value);

    Bool is_single() const override;
    Bool is_metric(ConfigurationPropertyPath const& path) const override;
    Bool is_configurable() const override;
    SizeType cardinality() const override;

    ConfigurationPropertyInterface* clone() const override;

    ConfigurationPropertyInterface* at(ConfigurationPropertyPath const& path) override;

    T const& get() const override;
    void set(T const& value) override;
    void set(List<T> const& values);
    void set_single(ConfigurationPropertyPath const& path, int integer_value) override;
    Map<ConfigurationPropertyPath,List<int>> integer_values() const override;
  protected:
    void local_set_single(int integer_value) override;
    List<int> local_integer_values() const override;
    List<SharedPointer<T>> values() const override;
  private:
    List<T> _values;
};

//! \brief A property that specifies a list of objects deriving from an interface \a T
//! \details T must define the clone() method to support interfaces.
template<class T> class InterfaceListConfigurationProperty final : public ConfigurationPropertyBase<T> {
  public:
    InterfaceListConfigurationProperty();
    InterfaceListConfigurationProperty(List<SharedPointer<T>> const& list);
    InterfaceListConfigurationProperty(T const& value);

    Bool is_single() const override;
    Bool is_metric(ConfigurationPropertyPath const& path) const override;
    Bool is_configurable() const override;
    SizeType cardinality() const override;

    ConfigurationPropertyInterface* clone() const override;

    ConfigurationPropertyInterface* at(ConfigurationPropertyPath const& path) override;

    T const& get() const override;
    Map<ConfigurationPropertyPath,List<int>> integer_values() const override;
    void set(T const& value) override;
    void set(SharedPointer<T> const& value);
    void set(List<SharedPointer<T>> const& values);
    void set_single(ConfigurationPropertyPath const& path, int integer_value) override;
  protected:
    void local_set_single(int integer_value) override;
    List<int> local_integer_values() const override;
    List<SharedPointer<T>> values() const override;
  private:
    List<SharedPointer<T>> _values;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_HPP
