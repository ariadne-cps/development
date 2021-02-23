/***************************************************************************
 *            set_wrapper.hpp
 *
 *  Copyright  2020  Pieter Collins
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

/*! \file set_wrapper.hpp
 *  \brief Wrappers for concrete sets
 */

#ifndef ARIADNE_SET_WRAPPER_HPP
#define ARIADNE_SET_WRAPPER_HPP

#include "utility/handle.hpp"
#include "set.hpp"

namespace Ariadne {

template<class SET, class P, class T> class OpenSetWrapper;
template<class SET, class P, class T> class ClosedSetWrapper;
template<class SET, class P, class T> class RegularSetWrapper;
template<class SET, class P, class T> class OvertSetWrapper;
template<class SET, class P, class T> class CompactSetWrapper;
template<class SET, class P, class T> class LocatedSetWrapper;
template<class SET, class P, class T> class RegularLocatedSetWrapper;

template<class SET, class P, class T> class LowerMeasurableSetWrapper;


template<class SET, class T> class OpenSetWrapper<SET,ValidatedTag,T>
    : public virtual ValidatedOpenSet<T>::Interface
    , public SET
{
    using P = ValidatedTag;
    SET const& base() const { return *this; }
  public:
    using typename OpenSetInterface<P,T>::DimensionType;
    using typename OpenSetInterface<P,T>::BasicSetType;

    using SET::SET;
    explicit OpenSetWrapper(SET const& set) : SET(set) { }
    virtual OutputStream& _write(OutputStream& os) const final override { return os << this->base(); }
    virtual OpenSetInterface<ValidatedTag,T>* clone() const final { return new OpenSetWrapper<SET,P,T>(this->base()); }
    virtual DimensionType dimension() const final override { return this->base().dimension(); }
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bs) const final { return this->base().overlaps(bs); }
    virtual ValidatedLowerKleenean covers(const BasicSetType& bs) const final { return this->base().covers(bs); }
    // Needed to prevent ambiguity
    friend OutputStream& operator<<(OutputStream& os, OpenSetWrapper<SET,P,T> ops) { return os << ops.base(); }
};

template<class SET, class T> class ClosedSetWrapper<SET,ValidatedTag,T>
    : public virtual ValidatedClosedSet<T>::Interface
    , public SET
{
    using P = ValidatedTag;
    SET const& base() const { return *this; }
  public:
    using typename ClosedSetInterface<P,T>::DimensionType;
    using typename ClosedSetInterface<P,T>::BasicSetType;

    using SET::SET;
    explicit ClosedSetWrapper(SET const& set) : SET(set) { }
    virtual OutputStream& _write(OutputStream& os) const final override { return os << this->base(); }
    virtual ClosedSetInterface<ValidatedTag,T>* clone() const final { return new ClosedSetWrapper<SET,P,T>(this->base()); }
    virtual DimensionType dimension() const final override { return this->base().dimension(); }
    virtual ValidatedLowerKleenean separated(const BasicSetType& bs) const final { return this->base().separated(bs); }
    // Needed to prevent ambiguity
    friend OutputStream& operator<<(OutputStream& os, ClosedSetWrapper<SET,P,T> ops) { return os << ops.base(); }
};

template<class SET, class T> class CompactSetWrapper<SET,ValidatedTag,T>
    : public virtual ValidatedCompactSet<T>::Interface
    , public SET
{
    using P = ValidatedTag;
    SET const& base() const { return *this; }
  public:
    using typename CompactSetInterface<P,T>::DimensionType;
    using typename CompactSetInterface<P,T>::BasicSetType;
    using typename CompactSetInterface<P,T>::BoundingSetType;

    using SET::SET;
    explicit CompactSetWrapper(SET const& set) : SET(set) { }
    virtual OutputStream& _write(OutputStream& os) const final override { return os << this->base(); }
    virtual CompactSetInterface<ValidatedTag,T>* clone() const final { return new CompactSetWrapper<SET,P,T>(this->base()); }
    virtual DimensionType dimension() const final override { return this->base().dimension(); }
#warning
    virtual BoundingSetType bounding_box() const final { return static_cast<BoundingSetType>(this->base().bounding_box()); }
    virtual ValidatedLowerKleenean inside(const BasicSetType& bs) const final { return this->base().inside(bs); }
    virtual ValidatedLowerKleenean separated(const BasicSetType& bs) const final { return this->base().separated(bs); }
    // Needed to prevent ambiguity
    friend OutputStream& operator<<(OutputStream& os, CompactSetWrapper<SET,P,T> ops) { return os << ops.base(); }
};

template<class SET, class T> class RegularSetWrapper<SET,ValidatedTag,T>
    : public virtual ValidatedRegularSet<T>::Interface
    , public SET
{
    using P = ValidatedTag;
    SET const& base() const { return *this; }
  public:
    using typename RegularSetInterface<P,T>::DimensionType;
    using typename RegularSetInterface<P,T>::BasicSetType;

    using SET::SET;
    explicit RegularSetWrapper(SET const& set) : SET(set) { }
    virtual OutputStream& _write(OutputStream& os) const final override { return os << this->base(); }
    virtual RegularSetInterface<ValidatedTag,T>* clone() const final { return new RegularSetWrapper<SET,P,T>(this->base()); }
    virtual DimensionType dimension() const final override { return this->base().dimension(); }
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bs) const final { return this->base().overlaps(bs); }
    virtual ValidatedLowerKleenean covers(const BasicSetType& bs) const final { return this->base().covers(bs); }
    virtual ValidatedLowerKleenean separated(const BasicSetType& bs) const final { return this->base().separated(bs); }
    // Needed to prevent ambiguity
    friend OutputStream& operator<<(OutputStream& os, RegularSetWrapper<SET,P,T> ops) { return os << ops.base(); }
};


template<class SET, class T> class RegularLocatedSetWrapper<SET,ValidatedTag,T>
    : public virtual RegularLocatedSet<ValidatedTag,T>::Interface
    , public SET
{
    using P = ValidatedTag;
    SET const& base() const { return *this; }
  public:
    using typename SetInterface<P,T>::DimensionType;
    using typename SetInterface<P,T>::BasicSetType;
    using typename SetInterface<P,T>::BoundingSetType;

    using SET::SET;
    RegularLocatedSetWrapper(SET const& set) : SET(set) { }
    virtual OutputStream& _write(OutputStream& os) const final override { return os << this->base(); }
    virtual RegularLocatedSetInterface<ValidatedTag,T>* clone() const final {
        return new RegularLocatedSetWrapper<SET,P,T>(this->base()); }
    virtual DimensionType dimension() const final override { return this->base().dimension(); }
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bs) const final { return this->base().overlaps(bs); }
    virtual ValidatedLowerKleenean covers(const BasicSetType& bs) const final { return this->base().covers(bs); }
    virtual ValidatedLowerKleenean separated(const BasicSetType& bs) const final { return this->base().separated(bs); }
    virtual ValidatedLowerKleenean inside(const BasicSetType& bs) const final { return this->base().inside(bs); }
    virtual BoundingSetType bounding_box() const final { return static_cast<BoundingSetType>(this->base().bounding_box()); }
    // Needed to prevent ambiguity
    friend OutputStream& operator<<(OutputStream& os, RegularLocatedSetWrapper<SET,P,T> ops) { return os << ops.base(); }
};

template<class SET> OpenSet<ValidatedTag,Real> wrap_open(SET const& s) {
    using P=ValidatedTag; using T=Real;
    return OpenSet<P,T>(std::make_shared<OpenSetWrapper<SET,P,T>>(s));
}

} // namespace Ariadne

#endif
