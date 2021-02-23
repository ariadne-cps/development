/***************************************************************************
 *            taylor_multifunction.hpp
 *
 *  Copyright 2012-21  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file taylor_multifunction.hpp
 *  \brief Multivalued functions based on Taylor polynomials
 */

#ifndef ARIADNE_TAYLOR_MULTIFUNCTION_HPP
#define ARIADNE_TAYLOR_MULTIFUNCTION_HPP


#include "taylor_model.hpp"
#include "scaled_function_patch.hpp"

#include "multifunction.hpp"


#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/algebra.hpp"
#include "../algebra/multi_index.hpp"
#include "../function/polynomial.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/evaluate.hpp"
#include "../geometry/interval.hpp"
#include "../geometry/box.hpp"
#include "../geometry/set_wrapper.hpp"


namespace Ariadne {

inline ValidatedLowerKleenean operator<(UpperBound<FloatMP> const& x1, LowerBound<FloatDP> const& x2) {
    if (x1.raw()<x2.raw()) { return true; } else { return ValidatedLowerKleenean(ValidatedKleenean(indeterminate)); }
}
inline ValidatedLowerKleenean operator>(LowerBound<FloatMP> const& x1, UpperBound<FloatDP> const& x2) {
    if (x1.raw()>x2.raw()) { return true; } else { return ValidatedLowerKleenean(ValidatedKleenean(indeterminate)); }
}


template<class M> class ScaledFunctionPatchMixin;

template<class U> class IntervalSetBase : public Interval<U> {
  protected:
    IntervalSetBase(Interval<U>&& ivl) : Interval<U>(ivl) { }
    template<class UU> decltype(auto) _inside(Interval<UU> const& ivl) const { return inside(*this,ivl); }
    template<class UU> decltype(auto) _separated(Interval<UU> const& ivl) const { return separated(*this,ivl); }
};

template<class U> class IntervalSet : public IntervalSetBase<U> {
  public:
    IntervalSet(Interval<U>&& ivl) : IntervalSetBase<U>(std::move(ivl)) { }
    Interval<U> const& bounding_box() const { return *this; }
    template<class UU> decltype(auto) inside(Interval<UU> const& ivl) const { return this->_inside(ivl); }
    template<class UU> decltype(auto) separated(Interval<UU> const& ivl) const { return this->_separated(ivl); }
};



template<class F> class ScaledFunctionPatchFactory<ValidatedIntervalTaylorModel<F>> {
  public:
    ScaledFunctionPatchFactory(Sweeper<F> const&) { }
};

template<class F> class ScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<F>>
    : public MultifunctionInterface<ValidatedTag, Real(RealVector),CompactSet>
{
    using M=ValidatedIntervalTaylorModel<F>; using P=ValidatedTag; using RES=Real;
    ScaledFunctionPatch<M> const& _upcast() const { return static_cast<ScaledFunctionPatch<M>const&>(*this); }
  protected:
    virtual CompactSet<P,RES> _call(Argument<Number<P>> const& x) const;
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_upcast(); }
};

template<class F> class VectorScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<F>>
    : public MultifunctionInterface<ValidatedTag, RealVector(RealVector), CompactSet>
{
    using M=ValidatedIntervalTaylorModel<F>; using P=ValidatedTag; using RES=RealVector;
    VectorScaledFunctionPatch<M> const& _upcast() const { return static_cast<VectorScaledFunctionPatch<M>const&>(*this); }
  protected:
    virtual CompactSet<P,RES> _call(Argument<Number<P>> const& x) const;
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_upcast(); }
};

template<class F> class ValidatedIntervalTaylorFunctionModel
    : public ScaledFunctionPatch<ValidatedIntervalTaylorModel<F>>
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<F>; using RES=Real;
  public:
    ValidatedIntervalTaylorFunctionModel(ScaledFunctionPatch<M> const& f) : ScaledFunctionPatch<M>(f) { }
    using ScaledFunctionPatch<M>::ScaledFunctionPatch;
    Interval<UpperBound<F>> operator() (Vector<Bounds<F>> const& x) const {
        return evaluate(this->model(),Ariadne::unscale(x,this->domain())); }
    CompactSet<P,RES> operator() (Vector<Number<P>> const& x) const { return this->_call(x); }
};

template<class F> class ValidatedVectorIntervalTaylorFunctionModel
    : public VectorScaledFunctionPatch<ValidatedIntervalTaylorModel<F>>
{
    using P=ValidatedTag; using M=ValidatedIntervalTaylorModel<F>; using RES=RealVector;
  public:
    ValidatedVectorIntervalTaylorFunctionModel(VectorScaledFunctionPatch<M> const& f) : VectorScaledFunctionPatch<M>(f) { }
    using VectorScaledFunctionPatch<M>::VectorScaledFunctionPatch;
    Box<Interval<UpperBound<F>>> operator() (Vector<Bounds<F>> const& x) const {
        return evaluate(this->models(),Ariadne::unscale(x,this->domain())); }
    CompactSet<P,RES> operator() (Vector<Number<P>> const& x) const { return this->_call(x); }
};


} // namespace Ariadne

#endif

