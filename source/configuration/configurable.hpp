/***************************************************************************
 *            configuration/configurable.hpp
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

/*! \file configuration/configurable.hpp
 *  \brief Classes for configurable tasks.
 *  \details If a class C needs a configuration, then it:
 *  1) must derive from Configurable<C>
 *  2) a Configuration<C> must be specified, which derives from SearchableConfiguration
 *  3) C constructors must take a Configuration<C> object
 *  For a base class B and one of its derivations D, D derives from B as usual, while the configuration design is a bit more complicated:
 *  1) B follows the design for a configurable class as described above
 *  2) Configuration<D> derives from Configuration<B>, with get/set accessors for Configuration<B> also, since Configuration<B> can't be passed
 *     (you could make Configuration<B> constructors protected, if you prefer)
 *  3) D needs its own configuration() override, which static casts B::configuration() onto Configuration<D>
 */

#ifndef ARIADNE_CONFIGURABLE_HPP
#define ARIADNE_CONFIGURABLE_HPP

#include "utility/pointer.hpp"

namespace Ariadne {

class SearchableConfiguration;

class ConfigurableInterface {
  public:
    virtual SearchableConfiguration const& searchable_configuration() const = 0;
};

//! \brief Base template class to be specialised while deriving from SearchableConfigurationInterface
template<class C> struct Configuration;

//! \brief Is-a component that provides a configuration
//! \details Since the configuration returned is const, a Configurable object should be constructed from
//! a pre-set configuration. If the configuration must specify certain properties or if some properties
//! must be coherent with the Configurable (e.g., the system used by a Configurable evolver), then a Builder
//! approach should be used for creation of the configuration, which should become an immutable object.
template<class C> class Configurable : public ConfigurableInterface {
    friend struct Configuration<C>;
  public:
    Configurable(Configuration<C> const& config);
    Configuration<C> const& configuration() const;
    SearchableConfiguration const& searchable_configuration() const override;
  private:
    SharedPointer<Configuration<C>> _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
