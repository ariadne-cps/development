/***************************************************************************
 *            concurrency/configurable.tpl.hpp
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

#ifndef ARIADNE_CONFIGURABLE_TPL_HPP
#define ARIADNE_CONFIGURABLE_TPL_HPP

#include "utility/macros.hpp"
#include "configurable.hpp"

namespace Ariadne {

template<class C> Configurable<C>::Configurable(Configuration<C> const& config) : _configuration(new Configuration<C>(config)) { }

template<class C> Configuration<C> const& Configurable<C>::configuration() const {
    return *_configuration;
}

template<class C> SearchableConfiguration const& Configurable<C>::searchable_configuration() const {
    return dynamic_cast<SearchableConfiguration const &>(*_configuration);
}

} // namespace Ariadne

#endif // ARIADNE_CONFIGURABLE_TPL_HPP
