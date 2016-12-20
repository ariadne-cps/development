/***************************************************************************
 *            function_model_mixin.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file function_model_mixin.h
 *  \brief Mixin for concrete functions on bounded domains.
 */

#ifndef ARIADNE_FUNCTION_MODEL_MIXIN_H
#define ARIADNE_FUNCTION_MODEL_MIXIN_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_model_interface.h"
#include "function/function_model.h"

#include "numeric/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/operations.h"
#include "geometry/box.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

namespace Ariadne {

template<class FM, class P, class PR=Precision64, class PRE=PR> class ScalarFunctionModelMixin;
template<class FM, class P, class PR=Precision64, class PRE=PR> class VectorFunctionModelMixin;

template<class FM, class P, class PR, class PRE> class ScalarFunctionModelMixin
    : public virtual ScalarFunctionModelInterface<P,PR,PRE>
    , public ScalarFunctionMixin<FM,P>
{
    typedef FloatError<PR> NormType;
  public:
    FM apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<P,PR,PRE>* _clone() const override {
        return new FM(static_cast<const FM&>(*this)); }
    NormType const _norm() const override {
        return norm(static_cast<const FM&>(*this)); }
    ScalarFunctionModelInterface<P,PR,PRE>* _antiderivative(SizeType j) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j)); }
    ScalarFunctionModelInterface<P,PR,PRE>* _antiderivative(SizeType j, CanonicalNumericType<P,PR,PRE> c) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j,c)); }
     ScalarFunctionModelInterface<P,PR,PRE>* _restriction(const ExactBoxType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    ScalarFunctionModelInterface<P,PR,PRE>* _apply(OperatorCode op) const override {
        return new FM(this->apply(op)); }
    CanonicalNumericType<P,PR,PRE> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterface<P,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
    ScalarFunctionModelInterface<P,PR,PRE>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const override {
        return new FM(embed(d1,static_cast<const FM&>(*this),d2)); }
    Boolean _refines(const ScalarFunctionModelInterface<P,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return refines(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    Boolean _inconsistent(const ScalarFunctionModelInterface<P,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return inconsistent(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    ScalarFunctionModelInterface<P,PR,PRE>* _refinement(const ScalarFunctionModelInterface<P,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return new FM(refinement(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Void _iadd(const CanonicalNumericType<P,PR,PRE>& c) override {
        static_cast<FM&>(*this)+=c; }
    Void _imul(const CanonicalNumericType<P,PR,PRE>& c) override {
        static_cast<FM&>(*this)*=c; }
    Void _isma(const CanonicalNumericType<P,PR,PRE>& c, const ScalarFunctionModelInterface<P,PR,PRE>& f) override {
        static_cast<FM&>(*this)+=c*dynamic_cast<const FM&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<P,PR,PRE>& f1, const ScalarFunctionModelInterface<P,PR,PRE>& f2) override {
        static_cast<FM&>(*this)+=dynamic_cast<const FM&>(f1)*dynamic_cast<const FM&>(f2); }
};

template<class FM, class P, class PR, class PRE> FM ScalarFunctionModelMixin<FM,P,PR,PRE>::apply(OperatorCode op) const {
    const FM& f=static_cast<const FM&>(*this);
    switch(op) {
        case OperatorCode::NEG: return neg(f);
        case OperatorCode::REC: return rec(f);
        case OperatorCode::EXP: return exp(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<P,PR,PRE>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}


template<class FM, class P, class PR, class PRE> class VectorFunctionModelMixin
    : public virtual VectorFunctionModelInterface<P,PR,PRE>
    , public VectorFunctionMixin<FM,P>
{
    typedef typename Element<FM>::Type ScalarFunctionType;
    typedef FloatError<PR> NormType;
  public:
    virtual VectorFunctionModelInterface<P,PR,PRE>* _clone() const override { return new FM(static_cast<const FM&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionModelInterface<P,PR,PRE>& sf) override {
        if(!dynamic_cast<const typename FM::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<FM&>(*this).FM::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    virtual VectorFunctionModelInterface<P,PR,PRE>* _derivative(SizeType j) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    NormType const _norm() const override {
         return norm(static_cast<const FM&>(*this)); }
    VectorFunctionModelInterface<P,PR,PRE>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const override {
        return heap_copy(embed(d1,static_cast<const FM&>(*this),d2)); }
    VectorFunctionModelInterface<P,PR,PRE>* _restriction(const ExactBoxType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    Void _adjoin(const ScalarFunctionModelInterface<P,PR,PRE>& f) override {
        static_cast<FM&>(*this).FM::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<P,PR,PRE>* _join(const VectorFunctionModelInterface<P,PR,PRE>& f) const override {
        return heap_copy(join(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    VectorFunctionModelInterface<P,PR,PRE>* _combine(const VectorFunctionModelInterface<P,PR,PRE>& f) const override {
        return heap_copy(combine(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Vector<CanonicalNumericType<P,PR,PRE>> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterface<P,PR,PRE>* _compose(const ScalarFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,PR,PRE>* _compose(const VectorFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    ScalarFunctionModelInterface<P,PR,PRE>* _unchecked_compose(const ScalarFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,PR,PRE>* _unchecked_compose(const VectorFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const FM&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
};

} // namespace Ariadne

#endif