/***************************************************************************
 *            function.h
 *
 *  Copyright 2008-12  Pieter Collins
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

/*! \file function.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/declarations.h"
#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"
#include "utility/metaprogramming.h"

#include "function/function_interface.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/differential.h"

#include "geometry/interval.h"
#include "geometry/box.h"

namespace Ariadne {

inline SizeOne dimension(IntervalDomain dom) { return SizeOne(); }
inline SizeType dimension(BoxDomain dom) { return dom.dimension(); }

template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;

template<class P, class D=BoxDomain> struct VectorFunctionElementReference;

template<class T> class Variable;
typedef Variable<Real> RealVariable;
template<class T> class Expression;
typedef Expression<Real> RealExpression;
template<class P, class D, class C> class FunctionExpression;

class RealDomain : public IntervalDomain {
  public:
    RealDomain() : IntervalDomain(-inf,+inf) { }
};

class EuclideanDomain : public BoxDomain {
  public:
    EuclideanDomain(SizeType n) : BoxDomain(n,RealDomain()) { }
};

template<class P>
class FunctionConstructors {
    static_assert(Or<IsSame<P,ApproximateTag>,IsSame<P,ValidatedTag>,IsSame<P,EffectiveTag>>::value,"P must be an information level/paradigm.");
    typedef Number<P> Y;
  public:
    typedef Y NumericType;

    static ScalarFunction<P,BoxDomain> zero(SizeType as);
    static ScalarFunction<P,BoxDomain> constant(SizeType as, NumericType c);
    static ScalarFunction<P,BoxDomain> coordinate(SizeType as, SizeType j);
    static List<ScalarFunction<P,BoxDomain>> coordinates(SizeType ns);
    static VectorFunction<P,BoxDomain> zeros(SizeType rs, SizeType as);
    static VectorFunction<P,BoxDomain> identity(SizeType ns);

    static ScalarFunction<P,BoxDomain> zero(BoxDomain dom);
    static ScalarFunction<P,BoxDomain> constant(BoxDomain dom, NumericType c);
    static ScalarFunction<P,BoxDomain> coordinate(BoxDomain dom, SizeType j);
    static List<ScalarFunction<P,BoxDomain>> coordinates(BoxDomain dom);
    static VectorFunction<P,BoxDomain> zeros(SizeType rs, BoxDomain dom);
    static VectorFunction<P,BoxDomain> identity(BoxDomain dom);

    static ScalarFunction<P,IntervalDomain> zero();
    static ScalarFunction<P,IntervalDomain> constant(NumericType c);
    static ScalarFunction<P,IntervalDomain> coordinate();
    static VectorFunction<P,IntervalDomain> zeros(SizeType rs);
    static ScalarFunction<P,IntervalDomain> identity();

    static ScalarFunction<P,IntervalDomain> zero(IntervalDomain dom);
    static ScalarFunction<P,IntervalDomain> constant(IntervalDomain dom, NumericType c);
    static ScalarFunction<P,IntervalDomain> coordinate(IntervalDomain dom, SizeType j);
    static VectorFunction<P,IntervalDomain> zeros(SizeType rs, IntervalDomain dom);
    static ScalarFunction<P,IntervalDomain> identity(IntervalDomain dom);

    static VectorFunction<P,BoxDomain> constant(BoxDomain dom, Vector<NumericType> c);
    static VectorFunction<P,IntervalDomain> constant(IntervalDomain dom, Vector<NumericType> c);

};

template<class P, class X> using EvaluateType = decltype(declval<ScalarFunctionInterface<P>>()._evaluate(declval<Vector<X>>()));

template<class P, class D, class C> class FunctionFacade {
};

template<class P> class FunctionFacade<P,IntervalDomain,IntervalDomain> {
  public:
    FunctionExpression<P,IntervalDomain,IntervalDomain> operator() (const RealVariable& x) const;
};

template<class P> class FunctionFacade<P,IntervalDomain,BoxDomain> {
  public:
    FunctionExpression<P,IntervalDomain,BoxDomain> operator() (const RealVariable& x) const;
};

template<class P> class FunctionFacade<P,BoxDomain,IntervalDomain> {
    typedef Number<P> Y;
  public:
    template<class X> Covector<EvaluateType<P,X>> gradient(Vector<X> const& x) const;
    FunctionExpression<P,BoxDomain,IntervalDomain> operator() (const Vector<RealVariable>& x) const;
    //FunctionExpression<P,BoxDomain,IntervalDomain> operator() (const Vector<RealExpression>& x) const;
};

template<class P> class FunctionFacade<P,BoxDomain,BoxDomain> {
    typedef Number<P> Y;
  public:
    template<class X> Matrix<EvaluateType<P,X>> jacobian(Vector<X> const& x) const;
    FunctionExpression<P,BoxDomain,BoxDomain> operator() (const Vector<RealVariable>& x) const;
};

template<class P, class D, class C> class DeclareFunctionOperations;
template<class P, class D> class DeclareFunctionOperations<P,D,IntervalDomain>
    : DeclareTranscendentalAlgebraOperations<Function<P,D,IntervalDomain>,Number<P>>
    , DeclareMixedArithmeticOperators<Function<P,D,IntervalDomain>,Int> { };
template<class P, class D> class DeclareFunctionOperations<P,D,BoxDomain>
    : DeclareVectorAlgebraOperators<Function<P,D,BoxDomain>,Function<P,D,IntervalDomain>,Number<P>> { };


//! \ingroup FunctionModule
//! \brief A generic scalar function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X\f$.
template<class P, class D, class C>
class Function
    : public FunctionConstructors<P>
    , public FunctionFacade<P,D,C>
    , public DeclareFunctionOperations<P,D,C>
{
    static_assert(IsStronger<P,ApproximateTag>::value,"P must be an information level/paradigm.");
    typedef Number<P> Y;
  protected:
    std::shared_ptr< const FunctionInterface<P,D,C> > _ptr;
  public:
    typedef P InformationTag;
    typedef P Paradigm;
    typedef Y NumericType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef decltype(dimension(declval<C>())) ResultSizeType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = typename ElementTraits<C>::template Type<Y>;

    explicit Function(EuclideanDomain dom);
    explicit Function(ResultSizeType rs, EuclideanDomain dom);

    explicit Function(DomainType dom);
    explicit Function(ResultSizeType rs, DomainType dom);

    explicit Function(DomainType dom, Result<Formula<Y>>const& e);
    explicit Function(ResultSizeType rs, ScalarFunction<P,D> sf);

    Function(InitializerList<ScalarFunction<P,D>> const& lsf);
    Function(List<ScalarFunction<P,D>> const& lsf);
    Function(Vector<ScalarFunction<P,D>> const& lsf);

    ScalarFunction<P,D> create_zero() const { return ScalarFunction<P,D>::zero(this->domain()); }
    ScalarFunction<P,D> create_constant(NumericType c) const { return ScalarFunction<P,D>::constant(this->domain(),c); }
    ScalarFunction<P,D> create_coordinate(SizeType j) const { return ScalarFunction<P,D>::coordinate(this->domain(),j); }
    VectorFunction<P,D> create_constant(Vector<NumericType> c) const { return VectorFunction<P,D>::constant(this->domain(),c); }

    Function();
    explicit Function(FunctionInterface<P,D,C>* p) : _ptr(p) { }
    explicit Function(SharedPointer<FunctionInterface<P,D,C>> p) : _ptr(p) { }
    Function(const FunctionInterface<P,D,C>& t) : _ptr(t._clone()) { }
    Function<P,D,C>& operator=(const FunctionInterface<P,D,C>& f) {
        _ptr=std::shared_ptr< FunctionInterface<P,D,C> >(f._clone()); return *this; }
    Function<P,D,C>& operator=(const Result<NumericType>& c) {
        return (*this)=this->create_constant(c); }

    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
    Function(const Function<PP,D,C>& f)
        : _ptr(std::dynamic_pointer_cast< const FunctionInterface<P,D,C> >(f.managed_pointer())) { }
    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
        Function<P,D,C>& operator=(Result<NumericType> const& c); // { return *this=this->create_constant(c); }

    shared_ptr< const FunctionInterface<P,D,C> > managed_pointer() const  { return _ptr; }
    const FunctionInterface<P,D,C>* raw_pointer() const  { return _ptr.operator->(); }
    const FunctionInterface<P,D,C>& reference() const  { return _ptr.operator*(); }
    operator const FunctionInterface<P,D,C>& () const { return _ptr.operator*(); }

    DomainType domain() const {
        return this->reference().domain(); }
    CodomainType codomain() const {
        return this->reference().codomain(); }
    SizeType argument_size() const {
        return this->reference().argument_size(); }
    SizeType result_size() const {
        return this->reference().result_size(); }

    template<class X> auto operator() (const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
    template<class X> auto evaluate(const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
    template<class X> friend auto evaluate(const Function<P,D,C>& f, const Argument<X>& x) -> decltype(f(x)) {
        return f(x); }

    friend VectorFunction<P,D> operator*(ScalarFunction<P,D> const&, Vector<Y> const&);

    Function<P,D,C> derivative(SizeType k) const {
        return Function<P,D,C>(this->reference()._derivative(k)); }
    friend Function<P,D,C> derivative(Function<P,D,C> const& f, SizeType k) {
        return Function<P,D,C>(f.reference()._derivative(k)); }

    template<class X> decltype(auto) differential(const Argument<X>& x, DegreeType d) const {
        return this->_ptr->_evaluate(Differential<EvaluateType<P,X>>::identity(d,x)); }

    Void set(SizeType i, ScalarFunction<P,D>);
    Function<P,D,IntervalDomain> get(SizeType i) const;
    const Function<P,D,IntervalDomain> operator[](SizeType i) const;
    VectorFunctionElementReference<P,D> operator[](SizeType i);

    friend OutputStream& operator<<(OutputStream& os, Function<P,D,C> const& f) { f._ptr->write(os); return os; }
};

template<class P, class D, class C> inline OutputStream&
operator<<(OutputStream& os, const Function<P,D,C>& f) {
    return f.write(os); }

template<class P, class D, class C, class X> inline decltype(auto)
evaluate(const Function<P,D,C>& f, const ElementType<D,X>& x) {
    return f(x); }

template<class P, class D, class C, class X> inline decltype(auto)
evaluate(const Function<P,D,C>& f, const Scalar<X>& x) {
    return f(x); }

template<class P, class D, class C, class X> inline decltype(auto)
evaluate(const Function<P,D,C>& f, const Vector<X>& x) {
    return f(x); }


template<class P, class C, class X> inline decltype(auto)
differential(const Function<P,IntervalDomain,C>& f, const X& x, DegreeType d) {
    return f.differential(x,d); }

template<class P, class C, class X> inline decltype(auto)
differential(const Function<P,BoxDomain,C>& f, const Vector<X>& x, DegreeType d) {
    return f.differential(x,d); }


template<class P, class X> Scalar<EvaluateType<P,X>>
derivative(const ScalarUnivariateFunction<P>& f, const Scalar<X>& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Vector<EvaluateType<P,X>>
tangent(const VectorUnivariateFunction<P>& f, const Scalar<X>& x) {
    return column(differential(f,x,1u).jacobian(),0u); }

template<class P, class X> Covector<EvaluateType<P,X>>
gradient(const ScalarMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Matrix<EvaluateType<P,X>>
jacobian(const VectorMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).jacobian(); }

template<class P> template<class X> Covector<EvaluateType<P,X>>
FunctionFacade<P,BoxDomain,IntervalDomain>::gradient(Vector<X> const& x) const {
    return Ariadne::gradient(static_cast<Function<P,BoxDomain,IntervalDomain>const&>(*this),x);
}

template<class P> template<class X> Matrix<EvaluateType<P,X>>
FunctionFacade<P,BoxDomain,BoxDomain>::jacobian(Vector<X> const& x) const {
    return Ariadne::jacobian(static_cast<Function<P,BoxDomain,BoxDomain>const&>(*this),x);
}


inline ValidatedScalarFunction& operator*=(ValidatedScalarFunction& sf, const ExactNumber& c) { return sf*=ValidatedNumber(c); }
EffectiveVectorFunction operator*(const EffectiveNumericType& c, const EffectiveVectorFunction& vf);

EffectiveScalarFunction embed(SizeType as1, const EffectiveScalarFunction& f2, SizeType as3);
EffectiveVectorFunction embed(SizeType as1, const EffectiveVectorFunction& f2, SizeType as3);

EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2);
EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveScalarFunction& f2);
EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);

EffectiveScalarFunction compose(const EffectiveScalarFunction& f, const EffectiveVectorFunction& g);
EffectiveVectorFunction compose(const EffectiveVectorFunction& f, const EffectiveVectorFunction& g);

EffectiveScalarFunction lie_derivative(const EffectiveScalarFunction& g, const EffectiveVectorFunction& f);

Formula<EffectiveNumericType> formula(const EffectiveScalarFunction& f);
Vector< Formula<EffectiveNumericType> > formula(const EffectiveVectorFunction& f);
//RealExpression evaluate(EffectiveScalarFunction const& f, Vector<RealVariable> const& vars);



ValidatedVectorFunction join(const ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2);
ValidatedVectorFunction join(const ValidatedScalarFunction& f1, const ValidatedVectorFunction& f2);
ValidatedVectorFunction join(const ValidatedVectorFunction& f1, const ValidatedScalarFunction& f2);
ValidatedVectorFunction join(const ValidatedVectorFunction& f1, const ValidatedVectorFunction& f2);
ValidatedScalarFunction compose(const ValidatedScalarFunction& f, const ValidatedVectorFunction& g);
ValidatedVectorFunction compose(const ValidatedVectorFunction& f, const ValidatedVectorFunction& g);



template<class P, class D>
struct VectorFunctionElementReference
    : DeclareFunctionOperations<P,D,IntervalDomain>
{
    typedef IntervalDomain SC; typedef BoxDomain VC;
    typedef VectorFunctionElementReference<P,D> SelfType;
    typedef Function<P,D,SC> ElementType;
    typedef typename ElementType::NumericType NumericType;
    Function<P,D,VC>& _vf; SizeType _i;
    VectorFunctionElementReference<P,D>(Function<P,D,VC>& vf, SizeType i) : _vf(vf), _i(i) { }
    template<class WP> operator Function<WP,D,SC> () const;
    Void operator=(const Function<P,D,SC>& sf);
    VectorFunctionElementReference<P,D>& operator=(const VectorFunctionElementReference<P,D>& sfr);
    D domain() const;
    template<class XX> XX evaluate(const Vector<XX> & x) const;
    template<class XX> XX operator()(const Vector<XX> & x) const;
};

template<class P, class D> template<class WP> inline
VectorFunctionElementReference<P,D>::operator Function<WP,D,SC> () const {
    return this->_vf.get(this->_i);
}

template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(ScalarFunction<P,D>& sf, SizeType i) {
    throw std::runtime_error(""); }
template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(VectorFunction<P,D>& vf, SizeType i) {
    return VectorFunctionElementReference<P,D>(vf,i); }

template<class P, class D, class C> inline VectorFunctionElementReference<P,D> Function<P,D,C>::operator[](SizeType i) {
    return make_element_reference(*this,i); }

template<class P, class D, class C> inline ScalarFunction<P,D> const Function<P,D,C>::operator[](SizeType i) const {
    return this->get(i); }

template<class P, class D> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionElementReference<P,D>& vfe) {
    return  os << static_cast< Function<P,D,IntervalDomain> >(vfe); }

template<class P, class D> inline Void VectorFunctionElementReference<P,D>::operator=(const Function<P,D,SC>& sf) {
    _vf.set(_i,sf); }
template<class P, class D> inline VectorFunctionElementReference<P,D>& VectorFunctionElementReference<P,D>::operator=(const VectorFunctionElementReference<P,D>& sfr) {
    _vf.set(_i,static_cast<Function<P,D,SC>>(sfr)); return *this; }
template<class P, class D> inline D VectorFunctionElementReference<P,D>::domain() const {
    return _vf.domain(); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::evaluate(const Vector<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::operator()(const Vector<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }



UpperIntervalType evaluate_range(ScalarFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x);
Vector<UpperIntervalType> evaluate_range(VectorFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x);
Vector<Differential<UpperIntervalType>> derivative_range(VectorFunction<ValidatedTag>const& f, const Vector<Differential<UpperIntervalType>>& x);
Covector<UpperIntervalType> gradient_range(ValidatedScalarFunction const& f, const Vector<UpperIntervalType>& x);
Matrix<UpperIntervalType> jacobian_range(ValidatedVectorFunction const& f, const Vector<UpperIntervalType>& x);



/*
inline Matrix<UpperIntervalType> jacobian(VectorFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian(VectorFunction<ValidatedTag>const& f, const Vector<ExactInterval>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian_range(VectorFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }

// FIXME: Needed to override templated gradient and jacobian
inline Covector<UpperIntervalType> gradient(ScalarFunction<EffectiveTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Covector<UpperIntervalType>>(gradient(f,reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Covector<UpperIntervalType> gradient(ScalarFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return gradient(f,static_cast<Vector<UpperIntervalType>>(x)); }
inline Matrix<UpperIntervalType> jacobian(VectorFunction<EffectiveTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian(VectorFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return jacobian(f,static_cast<Vector<UpperIntervalType>>(x)); }
*/

template<class P> class FunctionFactory;
typedef FunctionFactory<ValidatedTag> ValidatedFunctionFactory;

template<>
class FunctionFactory<ValidatedTag>
{
    std::shared_ptr< const FunctionFactoryInterface<ValidatedTag> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>* ptr) : _ptr(ptr) { }
    FunctionFactory(std::shared_ptr< const FunctionFactoryInterface<ValidatedTag> > ptr) : _ptr(ptr) { }
    inline ValidatedScalarFunction create(const ExactBoxType& d, const ValidatedScalarFunctionInterface& f) const;
    inline ValidatedVectorFunction create(const ExactBoxType& d, const ValidatedVectorFunctionInterface& f) const;
    inline ValidatedScalarFunction create_zero(const ExactBoxType& d) const;
    inline ValidatedVectorFunction create_identity(const ExactBoxType& d) const;
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<ValidatedTag>& factory);
};

inline ValidatedScalarFunction FunctionFactoryInterface<ValidatedTag>::create(const ExactBoxType& domain, const ValidatedScalarFunctionInterface& function) const {
    return ValidatedScalarFunction(SharedPointer<ValidatedScalarFunctionInterface>(this->_create(domain,function))); }
inline ValidatedVectorFunction FunctionFactoryInterface<ValidatedTag>::create(const ExactBoxType& domain, const ValidatedVectorFunctionInterface& function) const {
    return ValidatedVectorFunction(SharedPointer<ValidatedVectorFunctionInterface>(this->_create(domain,function))); }
inline ValidatedScalarFunction FunctionFactoryInterface<ValidatedTag>::create_zero(const ExactBoxType& domain) const {
    return this->create(domain,EffectiveScalarFunction::zero(domain.dimension())); }
inline ValidatedVectorFunction FunctionFactoryInterface<ValidatedTag>::create_identity(const ExactBoxType& domain) const {
    return this->create(domain,EffectiveVectorFunction::identity(domain.dimension())); }

inline ValidatedScalarFunction FunctionFactory<ValidatedTag>::create(const ExactBoxType& domain, const ValidatedScalarFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedVectorFunction FunctionFactory<ValidatedTag>::create(const ExactBoxType& domain, const ValidatedVectorFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedScalarFunction FunctionFactory<ValidatedTag>::create_zero(const ExactBoxType& domain) const {
    return this->_ptr->create_zero(domain); }
inline ValidatedVectorFunction FunctionFactory<ValidatedTag>::create_identity(const ExactBoxType& domain) const {
    return this->_ptr->create_identity(domain); }

inline OutputStream& operator<<(OutputStream& os, const ValidatedFunctionFactory& factory) {
    factory._ptr->write(os); return os; }

} // namespace Ariadne

#endif
