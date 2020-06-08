/***************************************************************************
 *            numeric/number_wrapper.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/number_wrapper.hpp
 *  \brief
 */



#ifndef ARIADNE_NUMBER_WRAPPER_HPP
#define ARIADNE_NUMBER_WRAPPER_HPP

#include <cassert>

#include "../utility/module.hpp"
#include "../numeric/paradigm.hpp"

#include "number_interface.hpp"

#include "number.hpp"

#include "logical.hpp"
#include "builtin.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float_value.hpp"
#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_upper_bound.hpp"
#include "float_lower_bound.hpp"
#include "float_approximation.hpp"
#include "float_error.hpp"

#include "../symbolic/templates.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

/************ Number *********************************************************/

// Declare fallbacks for Integer, Dyadic and Rational numbers
Real sqrt(Real const&);
Real exp(Real const&);
Real log(Real const&);
Real sin(Real const&);
Real cos(Real const&);
Real tan(Real const&);
Real asin(Real const&);
Real acos(Real const&);
Real atan(Real const&);

// Declare fallbacks for use by Lower/UpperBound
Approximation<FloatDP> div(Approximation<FloatDP> const& x1, Approximation<FloatDP> const& x2);
Approximation<FloatDP> sqr(Approximation<FloatDP> const& x);
Approximation<FloatDP> rec(Approximation<FloatDP> const& x);
Approximation<FloatDP> sin(Approximation<FloatDP> const& x);
Approximation<FloatDP> cos(Approximation<FloatDP> const& x);
Approximation<FloatDP> tan(Approximation<FloatDP> const& x);
Approximation<FloatDP> asin(Approximation<FloatDP> const& x);
Approximation<FloatDP> acos(Approximation<FloatDP> const& x);

Approximation<FloatMP> div(Approximation<FloatMP> const& x1, Approximation<FloatMP> const& x2);
Approximation<FloatMP> sqr(Approximation<FloatMP> const& x);
Approximation<FloatMP> rec(Approximation<FloatMP> const& x);
Approximation<FloatMP> sin(Approximation<FloatMP> const& x);
Approximation<FloatMP> cos(Approximation<FloatMP> const& x);
Approximation<FloatMP> tan(Approximation<FloatMP> const& x);
Approximation<FloatMP> asin(Approximation<FloatMP> const& x);
Approximation<FloatMP> acos(Approximation<FloatMP> const& x);


template<class X> class NumberMixin;
template<class X> class NumberWrapper;

template<class X> inline X const* extract(NumberInterface const* y) {
     return dynamic_cast<NumberWrapper<X>const*>(y);
}

inline OutputStream& operator<<(OutputStream& os, NumberInterface const& y) { return y._write(os); }


// FIXME: Should test for other potential infinities
inline Comparison cmp(NumberInterface const& y1, NumberInterface const& y2) {
    Comparison res;
    FloatDPValue const* x1=extract<FloatDPValue>(&y1);
    FloatDPValue const* x2=extract<FloatDPValue>(&y2);
    if(x1) {
        if(x2) { res= cmp(ExactDouble(x1->raw().get_d()),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(ExactDouble(x1->raw().get_d()),y2._get_q()); }
    } else {
        if(x2) { res= cmp(y1._get_q(),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(y1._get_q(),y2._get_q()); }
    }
    return res;
}

template<class... YS> struct Aware;

template<class DI> struct Managed;
template<class DI> using ManagedTypes = typename Managed<DI>::Types;

template<class T> struct DispatcherTraits;
template<class T> using DispatcherInterface = typename DispatcherTraits<T>::Interface;



template<class I, class DI, class OP, class YS=ManagedTypes<DI>> struct RightOperableInterface;
template<class I, class DI, class OP, class Y> struct RightOperableInterface<I,DI,OP,Aware<Y>> : public virtual I {
    virtual I* _concrete_apply_right(OP op, Y const& y) const = 0;
};
template<class I, class DI, class OP, class Y, class... YS> struct RightOperableInterface<I,DI,OP,Aware<Y,YS...>>
    : public virtual RightOperableInterface<I,DI,OP,Aware<YS...>>
{
    using RightOperableInterface<I,DI,OP,Aware<YS...>>::_concrete_apply_right;
    virtual I* _concrete_apply_right(OP op, Y const& y) const = 0;
};

template<class I, class DI, class OP, class YS=ManagedTypes<DI>> struct OperableInterface;
template<class I, class DI, class OP, class Y> struct OperableInterface<I,DI,OP,Aware<Y>>
    : public virtual RightOperableInterface<I,DI,OP>
{
    virtual I* _concrete_apply_left(OP op, Y const& y) const = 0;
};
template<class I, class DI, class OP, class Y, class... YS> struct OperableInterface<I,DI,OP,Aware<Y,YS...>>
    : public virtual OperableInterface<I,DI,OP,Aware<YS...>>
{
    using OperableInterface<I,DI,OP,Aware<YS...>>::_concrete_apply_left;
    virtual I* _concrete_apply_left(OP op, Y const& y) const = 0;
};


/*
template<class OP, class X1, class X2> struct CanApply {
    template<class OPP, class XX1, class XX2, class=decltype(declval<OPP>()(declval<XX1>(),declval<XX2>()))>
        static std::true_type test(int);
    template<class OPP, class XX1, class XX2>
        static std::false_type test(...);
    static const bool value = decltype(test<OP,X1,X2>(1))::value;
};
template<class OP, class X1, class X2> NumberInterface* _concrete_apply(OP op, X1 const& x1, X2 const& x2) {
    if constexpr (CanApply<OP,X1,X2>::value) {
        return _make_number_wrapper(op(x1,x2));
    } else {
        String yc1=class_name<X1>(); String yc2=class_name<X2>();
        ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<x1<<", y2="<<x2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
    }
}
*/

template<class R> NumberInterface* _make_number_wrapper(R const& r) { return new NumberWrapper<R>(r); }

template<class OP, class X1, class X2> inline NumberInterface* _concrete_apply(OP op, X1 const& x1, X2 const& x2) {
    return _make_number_wrapper(op(x1,x2));
}
// Since OP(Value<F>,Value<F>) usually returns Bounds<F>, which is Validated, while an Exact answer is expected,
// disable binary operation on value
// FIXME: Prefer symbolic dispatch
template<class OP, class F> inline NumberInterface* _concrete_apply(OP op, Value<F> const& x1, Value<F> const& x2) {
    String yc1=class_name<Value<F>>(); String yc2=class_name<Value<F>>();
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<x1<<", y2="<<x2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}
template<class OP, class F> inline NumberInterface* _concrete_apply(OP op, Value<F> const& x1, Rational const& q2) {
    String yc1=class_name<Value<F>>(); String yc2=class_name<Rational>();
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<x1<<", y2="<<q2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}
template<class OP, class F> inline NumberInterface* _concrete_apply(OP op, Rational const& q1, Value<F> const& x2) {
    String yc1=class_name<Rational>(); String yc2=class_name<Value<F>>();
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<q1<<", y2="<<x2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}


template<class X, class I, class DI, class OP, class Y=ManagedTypes<DI>> struct RightOperableMixin;
template<class X, class I, class DI, class OP> struct RightOperableMixin<X,I,DI,OP,Aware<>>
    : virtual RightOperableInterface<I,DI,OP>
{
    X const& self() const { return static_cast<NumberWrapper<X>const&>(*this); }
};
template<class X, class I, class DI, class OP, class Y, class... YS> struct RightOperableMixin<X,I,DI,OP,Aware<Y,YS...>>
    : RightOperableMixin<X,I,DI,OP,Aware<YS...>>
{
    using RightOperableMixin<X,I,DI,OP,Aware<YS...>>::_concrete_apply_right;
    virtual I* _concrete_apply_right(OP op, Y const& y) const { return _concrete_apply(op,y,this->self()); }
};

template<class X, class I, class DI, class OP, class Y=ManagedTypes<DI>> struct OperableMixin;
template<class X, class I, class DI, class OP> struct OperableMixin<X,I,DI,OP,Aware<>>
    : virtual OperableInterface<I,DI,OP>
{
    X const& self() const { return static_cast<NumberWrapper<X>const&>(*this); }
};
template<class X, class I, class DI, class OP, class Y, class... YS> struct OperableMixin<X,I,DI,OP,Aware<Y,YS...>>
    : OperableMixin<X,I,DI,OP,Aware<YS...>>
{
    using OperableMixin<X,I,DI,OP,Aware<YS...>>::_concrete_apply_left;
    using OperableMixin<X,I,DI,OP,Aware<YS...>>::_concrete_apply_right;
    virtual I* _concrete_apply_left(OP op, Y const& y) const { return _concrete_apply(op,this->self(),y); }
    virtual I* _concrete_apply_right(OP op, Y const& y) const { return _concrete_apply(op,y,this->self()); }
};




class AlgebraicNumberInterface;
template<class F> class ConcreteNumberInterface;
template<class F, class FE> class ConcreteBallInterface;

template<> struct Managed<AlgebraicNumberInterface> {
    typedef Aware<Integer,Dyadic,Rational,Real> Types;
};
template<class F> struct Managed<ConcreteNumberInterface<F>> {
    typedef Aware<Value<F>,Ball<F>,Bounds<F>,UpperBound<F>,LowerBound<F>,Approximation<F>> Types;
};
template<class F, class FE> struct Managed<ConcreteBallInterface<F,FE>> {
    typedef Aware<Ball<F,FE>> Types;
};

template<> struct DispatcherTraits<Integer> { typedef AlgebraicNumberInterface Interface; };
template<> struct DispatcherTraits<Dyadic> { typedef AlgebraicNumberInterface Interface; };
template<> struct DispatcherTraits<Rational> { typedef AlgebraicNumberInterface Interface; };
template<> struct DispatcherTraits<Real> { typedef AlgebraicNumberInterface Interface; };

template<class F> struct DispatcherTraits<Value<F>> { typedef ConcreteNumberInterface<F> Interface; };
template<class F> struct DispatcherTraits<Ball<F>> { typedef ConcreteNumberInterface<F> Interface; };
template<class F> struct DispatcherTraits<Bounds<F>> { typedef ConcreteNumberInterface<F> Interface; };
template<class F> struct DispatcherTraits<UpperBound<F>> { typedef ConcreteNumberInterface<F> Interface; };
template<class F> struct DispatcherTraits<LowerBound<F>> { typedef ConcreteNumberInterface<F> Interface; };
template<class F> struct DispatcherTraits<Approximation<F>> { typedef ConcreteNumberInterface<F> Interface; };

template<class F, class FE> struct DispatcherTraits<Ball<F,FE>> { typedef ConcreteBallInterface<F,FE> Interface; };


class AlgebraicNumberInterface
    : public virtual RightOperableInterface<NumberInterface,AlgebraicNumberInterface,BinaryElementaryOperator>
{
};

template<class F> class ConcreteNumberInterface
    : public virtual RightOperableInterface<NumberInterface,ConcreteNumberInterface<F>,BinaryElementaryOperator>
    , public virtual OperableInterface<NumberInterface,AlgebraicNumberInterface,BinaryElementaryOperator>
{
};

template<class F, class FE> class ConcreteBallInterface
    : public virtual RightOperableInterface<NumberInterface,ConcreteBallInterface<F,FE>,BinaryElementaryOperator>
    , public virtual OperableInterface<NumberInterface,AlgebraicNumberInterface,BinaryElementaryOperator>
//    , public virtual OperableInterface<NumberInterface,ConcreteNumberInterface<F>,BinaryElementaryOperator>
{
};



template<class X, class DI=DispatcherInterface<X>> class ElementaryBinaryNumberDispatcherMixin;

template<class Y> class ElementaryBinaryNumberDispatcherMixin<Y,AlgebraicNumberInterface>
    : public RightOperableMixin<Y,NumberInterface,AlgebraicNumberInterface,BinaryElementaryOperator>
{
};

template<class Y, class F> class ElementaryBinaryNumberDispatcherMixin<Y,ConcreteNumberInterface<F>>
    : public RightOperableMixin<Y,NumberInterface,ConcreteNumberInterface<F>,BinaryElementaryOperator>
    , public OperableMixin<Y,NumberInterface,AlgebraicNumberInterface,BinaryElementaryOperator>
{
};

template<class Y, class F, class FE> class ElementaryBinaryNumberDispatcherMixin<Y,ConcreteBallInterface<F,FE>>
    : public RightOperableMixin<Y,NumberInterface,ConcreteBallInterface<F,FE>,BinaryElementaryOperator>
    , public OperableMixin<Y,NumberInterface,AlgebraicNumberInterface,BinaryElementaryOperator>
//    , public OperableMixin<Y,NumberInterface,ConcreteNumberInterface<F>,BinaryElementaryOperator>
{
};









template<class I, class X, class OP> inline I* _apply(X const& self, OP op, I const* self_ptr, I const* other_ptr) {
    using DI=DispatcherInterface<X>;
    auto aware_other_ptr=dynamic_cast<RightOperableInterface<I,DI,OP> const*>(other_ptr);
    if(aware_other_ptr) { return aware_other_ptr->_concrete_apply_right(op,self); }
    else { return other_ptr->_rapply(op,self_ptr); }
}
template<class I, class X, class OP> inline I* _rapply(X const& self, OP op, I const* self_ptr, I const* other_ptr) {
    using DI=DispatcherInterface<X>;
    auto aware_other_ptr=dynamic_cast<OperableInterface<I,DI,OP> const*>(other_ptr);
    if(aware_other_ptr) { return aware_other_ptr->_concrete_apply_left(op,self); }
    else { return make_symbolic(op,other_ptr,self_ptr); }
}







template<class OP> inline NumberInterface* make_symbolic(OP op, NumberInterface const* yp1, NumberInterface const* yp2) {
//    Handle<NumberInterface> y1(const_cast<NumberInterface*>(yp1)->shared_from_this());
//    Handle<NumberInterface> y2(const_cast<NumberInterface*>(yp2)->shared_from_this());
    String yc1=yp1->_class_name(); String yc2=yp2->_class_name();
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<*yp1<<", y2="<<*yp2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}




/*
template<class I, class OP> struct UnaryOperationInterface {
    virtual ~UnaryOperationInterface() = default;
    virtual I* _apply(OP op) const = 0;
};

template<class I, class OP> struct BinaryOperationInterface {
    virtual ~BinaryOperationInterface() = default;
    virtual I* _apply(OP op, I const* other) const = 0;
    virtual I* _rapply(OP op, I const* other) const = 0;
};

template<class I, class OP, class N> struct GradedOperationInterface {
    virtual ~GradedOperationInterface() = default;
    virtual I* _apply(OP op, N n) const = 0;
};

*/

template<class X, class I, class OP> struct UnaryOperationMixin : public virtual I {
    static X const& _cast(UnaryOperationMixin<X,I,OP> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(OP op) const final { return _make_wrapper(op(_cast(*this))); }
};

template<class X, class I, class OP> struct BinaryOperationMixin : public virtual I {
    static X const& _cast(BinaryOperationMixin<X,I,OP> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(OP op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(OP op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};

template<class X, class I, class OP, class N> struct GradedOperationMixin : public virtual I {
    static X const& _cast(GradedOperationMixin<X,I,OP,N> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(OP op, N n) const final { return _make_wrapper(op(_cast(*this),n)); }
};

template<class X> struct ElementaryUnaryNumberOperationsMixin
    : public UnaryOperationMixin<X,NumberInterface,UnaryElementaryOperator> { };
template<class X> struct ElementaryBinaryNumberOperationsMixin
    : public BinaryOperationMixin<X,NumberInterface,BinaryElementaryOperator> { };
template<class X,class N=Int> struct ElementaryGradedNumberOperationsMixin
    : public GradedOperationMixin<X,NumberInterface,GradedElementaryOperator,N> { };


template<class X> class NumberGetterMixin : public virtual NumberInterface {
  public:
  //  operator X const& () const { return static_cast<NumberMixin<X>const&>(*this); }
    static X const& _cast(NumberGetterMixin<X> const& self) { return static_cast<NumberWrapper<X>const&>(static_cast<NumberMixin<X>const&>(self)); }
    static X& _cast(NumberGetterMixin<X>& self) { return static_cast<NumberWrapper<X>&>(static_cast<NumberMixin<X>&>(self)); }

    typedef Paradigm<X> P;
    friend class Number<P>;

    virtual NumberInterface* _copy() const override { return new NumberWrapper<X>(_cast(*this)); }
    virtual NumberInterface* _move() override { return new NumberWrapper<X>(std::move(_cast(*this))); }

    // FIXME: Proper comparisons for ExactNumber.
    virtual LogicalValue _equals(NumberInterface const& y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y._paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,y)==Comparison::EQUAL ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(OrderTag(),dp) == y._get(OrderTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) == y._get(ApproximateTag(),dp)); }
    }
    virtual LogicalValue _less(NumberInterface const& y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y._paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,y)==Comparison::LESS ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        else if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(OrderTag(),dp) < y._get(OrderTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) < y._get(ApproximateTag(),dp));
        }
    }

    virtual Rational _get_q() const override {
        return this->_get_as<Rational>(); }

    virtual FloatDPBall _get(MetricTag,DoublePrecision pr,DoublePrecision pre) const override {
        return this->_get_as<FloatDPBall>(pr); }
    virtual FloatDPBounds _get(OrderTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPBounds>(pr); }
    virtual FloatDPUpperBound _get(UpperTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPUpperBound>(pr); }
    virtual FloatDPLowerBound _get(LowerTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPLowerBound>(pr); }
    virtual FloatDPApproximation _get(ApproximateTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPApproximation>(pr); }
    virtual FloatMPDPBall _get(MetricTag,MultiplePrecision pr, DoublePrecision pre) const override {
        return this->_get_as<FloatMPDPBall>(pr,pre); }
    virtual FloatMPBall _get(MetricTag, MultiplePrecision pr, MultiplePrecision pre) const override {
        return this->_get_as<FloatMPBall>(pr,pre); }
    virtual FloatMPBounds _get(OrderTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPBounds>(pr); }
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPUpperBound>(pr); }
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPLowerBound>(pr); }
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPApproximation>(pr); }

    virtual ParadigmCode _paradigm() const override { return P::code(); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R> inline R _get_as() const {
        if constexpr (Constructible<R,X>) { return static_cast<R>(_cast(*this)); }
        else { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << "\n"; throw ParadigmError(); } }
    template<class R, class PR> inline R _get_as(PR pr) const {
        if constexpr (Constructible<R,X,PR>) { return R(_cast(*this),pr); }
        else { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision " << pr << "\n"; throw ParadigmError(); } }
    template<class R, class PR, class PRE> inline R _get_as(PR pr, PRE pre) const {
        if constexpr (Constructible<R,X,PR,PRE>) { return R(_cast(*this),pr,pre); }
        else { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision " << pr << " and error precision " << pre << "\n"; throw ParadigmError(); } }
};

template<class X> class NumberMixin
    : public ElementaryBinaryNumberOperationsMixin<X>
    , public ElementaryUnaryNumberOperationsMixin<X>
    , public ElementaryGradedNumberOperationsMixin<X>
    , public ElementaryBinaryNumberDispatcherMixin<X>
    , public NumberGetterMixin<X>
{
  public:
    operator X const& () const { return static_cast<NumberWrapper<X>const&>(*this); }
    operator X& () { return static_cast<NumberWrapper<X>&>(*this); }
};


template<class X> class NumberWrapper
    : public X, public NumberMixin<X>
{
    inline static const X& _cast(const NumberWrapper<X>& x) { return static_cast<const X&>(x); }
    static_assert(not Same<X,Handle<NumberInterface>>,"X must be a concrete number, not a handle");
    static_assert(not Same<X,Number<Paradigm<X>>>,"X must be a concrete number, not a generic number");
  public:
    NumberWrapper(const X& a) : X(a) { }
    NumberWrapper(X&& a) : X(std::forward<X>(a)) { }
};


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_HPP */
