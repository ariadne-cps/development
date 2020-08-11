/***************************************************************************
 *            function/domain.hpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

/*! \file function/domain.hpp
 *  \brief Interval and box domains for functions.
 */

#ifndef ARIADNE_DOMAIN_HPP
#define ARIADNE_DOMAIN_HPP

#include "geometry/interval.hpp"
#include "geometry/box.hpp"

namespace Ariadne {

using IntervalDomainType = Interval<FloatDPValue>;
using BoxDomainType = Box<Interval<FloatDPValue>>;

class RealDomain {
  public:
    typedef SizeOne DimensionType;

    constexpr RealDomain() { }
    constexpr RealDomain(SizeOne) { }
    constexpr SizeOne dimension() const { return SizeOne(); }
    operator IntervalDomainType() const { return IntervalDomainType(-inf,+inf); }
    friend RealDomain intersection(RealDomain const& dom1, RealDomain const& dom2) { return RealDomain(); }
    friend Bool operator==(RealDomain const& dom1, RealDomain const& dom2) { return true; }
    friend Bool operator==(RealDomain const& dom1, IntervalDomainType const& dom2) { return IntervalDomainType(dom1)==dom2; }
    friend Bool operator==(IntervalDomainType const& dom1, RealDomain const& dom2) { return dom1==IntervalDomainType(dom2); }
    friend OutputStream& operator<<(OutputStream& os, RealDomain const& dom) { return os << "R"; }
};

class EuclideanDomain {
    SizeType _dim;
  public:
    typedef SizeType DimensionType;
    constexpr EuclideanDomain(SizeType dim) : _dim(dim) { }
    constexpr EuclideanDomain(SizeType dim, RealDomain) : _dim(dim) { }
    constexpr SizeType dimension() const { return this->_dim; }
    constexpr RealDomain operator[](SizeType ind) { return RealDomain(); }
    operator BoxDomainType() const { return BoxDomainType(this->dimension(),IntervalDomainType(RealDomain())); }
    friend EuclideanDomain intersection(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { assert(dom1==dom2); return dom1; }
    friend EuclideanDomain product(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { return EuclideanDomain(dom1.dimension()+dom2.dimension()); }
    friend EuclideanDomain product(EuclideanDomain const& dom1, RealDomain const& dom2) { return EuclideanDomain(dom1.dimension()+dom2.dimension()); }
    friend Bool operator==(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { return dom1.dimension() == dom2.dimension(); }
    friend Bool operator==(EuclideanDomain const& dom1, BoxDomainType const& dom2) { return BoxDomainType(dom1)==dom2; }
    friend Bool operator==(BoxDomainType const& dom1, EuclideanDomain const& dom2) { return dom1==BoxDomainType(dom2); }

    friend OutputStream& operator<<(OutputStream& os, EuclideanDomain const& dom) { return os << "R" << dom.dimension(); }
};


class UnitInterval {
  public:
    constexpr UnitInterval() { }
    constexpr SizeOne dimension() const { return SizeOne(); }
    operator IntervalDomainType() const { return IntervalDomainType(-1,+1); }
    friend UnitInterval intersection(UnitInterval const& dom1, UnitInterval const& dom2) { return UnitInterval(); }
    friend Bool operator==(UnitInterval const& dom1, UnitInterval const& dom2) { return true; }
    friend OutputStream& operator<<(OutputStream& os, UnitInterval const& dom) { return os << "[-1:+1]"; }
};

class UnitBox {
    SizeType _dim;
  public:
    constexpr UnitBox(SizeType dim) : _dim(dim) { }
    constexpr UnitBox(SizeType dim, UnitInterval) : _dim(dim) { }
    constexpr SizeType dimension() const { return this->_dim; }
    constexpr UnitInterval operator[](SizeType ind) { return UnitInterval(); }
    operator BoxDomainType() const { return BoxDomainType(this->dimension(),IntervalDomainType(UnitInterval())); }
    friend UnitBox intersection(UnitBox const& dom1, UnitBox const& dom2) { assert(dom1==dom2); return dom1; }
    friend UnitBox product(UnitBox const& dom1, UnitBox const& dom2) { return UnitBox(dom1.dimension()+dom2.dimension()); }
    friend UnitBox product(UnitBox const& dom1, UnitInterval const& dom2) { return UnitBox(dom1.dimension()+dom2.dimension()); }
    friend Bool operator==(UnitBox const& dom1, UnitBox const& dom2) { return dom1.dimension() == dom2.dimension(); }
    friend OutputStream& operator<<(OutputStream& os, UnitBox const& dom) { return os << "[-1:+1]^" << dom.dimension(); }
};

struct ScalarTraits {
    template<class X> using Type=Scalar<X>;
    typedef Scalar<Real> Kind;
    typedef SizeOne SizeType;
    typedef IndexZero IndexType;
    using EntireDomainType = RealDomain;
    using BoundedDomainType = IntervalDomainType;
    template<class PR> using RangeType = Interval<FloatUpperBound<PR>>;
};

struct VectorTraits {
    template<class X> using Type=Vector<X>;
    typedef Vector<Real> Kind;
    typedef Ariadne::SizeType SizeType;
    typedef Ariadne::SizeType IndexType;
    using EntireDomainType = EuclideanDomain;
    using BoundedDomainType = BoxDomainType;
    template<class PR> using RangeType = Box<Interval<FloatUpperBound<PR>>>;
};

template<class R> struct DomainTraits;
template<> struct DomainTraits<RealScalar> : ScalarTraits { };
template<> struct DomainTraits<RealVector> : VectorTraits { };

template<class S> struct ElementTraits;
template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;

template<class UB> struct ElementTraits<Interval<UB>> : ScalarTraits { };
template<class IVL> struct ElementTraits<Box<IVL>> : VectorTraits { };
template<> struct ElementTraits<RealDomain> : ScalarTraits { };
template<> struct ElementTraits<EuclideanDomain> : VectorTraits { };
template<> struct ElementTraits<UnitInterval> : ScalarTraits { };
template<> struct ElementTraits<UnitBox> : VectorTraits { };

template<class... TS> using CartesianProductType = decltype(product(declval<TS>()...));

} // namespace Ariadne


#endif /* ARIADNE_DOMAIN_HPP */
