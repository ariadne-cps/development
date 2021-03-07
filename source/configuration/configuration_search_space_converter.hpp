/***************************************************************************
 *            configuration/search_space_converter.hpp
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

/*! \file configuration/configuration_search_space_converter.hpp
 *  \brief Classes for converting configuration property data from/to the integer search space.
 */

#ifndef ARIADNE_CONFIGURATION_SEARCH_SPACE_CONVERTER_HPP
#define ARIADNE_CONFIGURATION_SEARCH_SPACE_CONVERTER_HPP

#include "numeric/builtin.hpp"
#include "numeric/real.hpp"

namespace Ariadne {

//! \brief Interface for conversion from/into the integer search space
template<class T> struct ConfigurationSearchSpaceConverterInterface {
    //! \brief Convert the \a value into an integer value
    virtual int to_int(T const& value) const = 0;
    //! \brief Convert the \a value into a double value
    virtual double to_double(T const& value) const = 0;
    //! \brief Convert from an integer value \a i into the original value
    virtual T from_int(int i) const = 0;
    //! \brief Convert from a double into the original value
    virtual T from_double(double value) const = 0;

    virtual ConfigurationSearchSpaceConverterInterface* clone() const = 0;
    virtual ~ConfigurationSearchSpaceConverterInterface() = default;
};

template<class T> struct SearchSpaceConverterBase : public ConfigurationSearchSpaceConverterInterface<T> {
  public:
    T from_int(int i) const override final { return this->from_double(i); }
};

template<class T> struct Log10SearchSpaceConverter;
template<class T> struct Log2SearchSpaceConverter;
template<class T> struct LinearSearchSpaceConverter;

template<> struct Log10SearchSpaceConverter<ExactDouble> : SearchSpaceConverterBase<ExactDouble> {
    int to_int(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<int>::max();
        else if (value == -inf) return std::numeric_limits<int>::min();
        else return std::round(log(value.get_d())/log(10.0)); }
    double to_double(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<double>::max();
        else if (value == -inf) return std::numeric_limits<double>::min();
        else return log(value.get_d())/log(10.0); }
    ExactDouble from_double(double value) const override { return ExactDouble(exp(log(10.0)*value)); }
    ConfigurationSearchSpaceConverterInterface* clone() const override { return new Log10SearchSpaceConverter(*this); }
};

template<> struct Log2SearchSpaceConverter<ExactDouble> : SearchSpaceConverterBase<ExactDouble> {
    int to_int(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<int>::max();
        else if (value == -inf) return std::numeric_limits<int>::min();
        return std::round(log(value.get_d())/log(2.0)); }
    double to_double(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<double>::max();
        else if (value == -inf) return std::numeric_limits<double>::min();
        return log(value.get_d())/log(2.0); }
    ExactDouble from_double(double value) const override { return ExactDouble(exp(log(2.0)*value)); }
    ConfigurationSearchSpaceConverterInterface* clone() const override { return new Log2SearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<ExactDouble> : SearchSpaceConverterBase<ExactDouble> {
    int to_int(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<int>::max();
        else if (value == -inf) return std::numeric_limits<int>::min();
        return round(value).get<int>(); }
    double to_double(ExactDouble const& value) const override { return value.get_d(); }
    ExactDouble from_double(double value) const override { return ExactDouble(value); }
    ConfigurationSearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<int> : SearchSpaceConverterBase<int> {
    int to_int(int const& value) const override { return value; }
    double to_double(int const& value) const override { return value; }
    int from_double(double i) const override { return std::round(i); }
    ConfigurationSearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<DegreeType> : SearchSpaceConverterBase<DegreeType> {
    int to_int(DegreeType const& value) const override { return value; }
    double to_double(DegreeType const& value) const override { return value; }
    DegreeType from_double(double i) const override { return DegreeType(i); }
    ConfigurationSearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_SEARCH_SPACE_CONVERTER_HPP
