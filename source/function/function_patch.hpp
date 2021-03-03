/***************************************************************************
 *            function/function_patch.hpp
 *
 *  Copyright  2011-20  Pieter Collins
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

/*! \file function/function_patch.hpp
 *  \brief Functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_PATCH_HPP
#define ARIADNE_FUNCTION_PATCH_HPP

#include "../numeric/number.decl.hpp"
#include "../function/function.decl.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"

#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/range.hpp"
#include "../algebra/operations.hpp"
#include "../algebra/algebra_interface.hpp"
#include "../algebra/algebra_mixin.hpp"

#include "../function/domain.hpp"
#include "../function/function_interface.hpp"
#include "../function/function_patch_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function.hpp"


namespace Ariadne {

template<class P, class D> struct AlgebraOperations<ScalarFunctionPatch<P,D>>;

// FIXME: Extend with univariate case
template<class P> class FunctionPatchFactory {
    SharedPointer<const FunctionPatchFactoryInterface<P>> _ptr;
    typedef RealScalar SARG;
    typedef RealVector VARG;
  public:
    typedef P Paradigm;
    typedef IntervalDomainType ScalarDomainType;
    typedef BoxDomainType VectorDomainType;

    operator const FunctionPatchFactoryInterface<P>& () const { return *_ptr; }

    explicit FunctionPatchFactory(const FunctionPatchFactoryInterface<P>* p) : _ptr(p) { }
    explicit FunctionPatchFactory(SharedPointer<const FunctionPatchFactoryInterface<P>> p) : _ptr(p) { }

    ScalarUnivariateFunctionPatch<P> create(ScalarDomainType const& dom, ScalarUnivariateFunction<P> const& f) const {
        return ScalarUnivariateFunctionPatch<P>(this->_ptr->_create(dom,f)); }
    VectorUnivariateFunctionPatch<P> create(ScalarDomainType const& dom, VectorUnivariateFunction<P> const& f) const {
        return VectorUnivariateFunctionPatch<P>(this->_ptr->_create(dom,f)); }
    ScalarMultivariateFunctionPatch<P> create(VectorDomainType const& dom, ScalarMultivariateFunction<P> const& f) const {
        return ScalarMultivariateFunctionPatch<P>(this->_ptr->_create(dom,f)); }
    VectorMultivariateFunctionPatch<P> create(VectorDomainType const& dom, VectorMultivariateFunction<P> const& f) const {
        return VectorMultivariateFunctionPatch<P>(this->_ptr->_create(dom,f)); }

    ScalarFunctionPatch<P,VARG> create_zero(VectorDomainType const& dom) const {
        return ScalarFunctionPatch<P,VARG>(this->_ptr->_create_zero(dom)); }
    ScalarFunctionPatch<P,VARG> create_constant(VectorDomainType const& dom, Number<P> const& c) const {
        return ScalarFunctionPatch<P,VARG>(this->_ptr->_create_constant(dom,c)); }
    ScalarFunctionPatch<P,VARG> create_coordinate(VectorDomainType const& dom, SizeType index) const {
        return ScalarFunctionPatch<P,VARG>(this->_ptr->_create_coordinate(dom,index)); }
    VectorFunctionPatch<P,VARG> create_zeros(SizeType rsize, VectorDomainType const& dom) const {
        return VectorFunctionPatch<P,VARG>(this->_ptr->_create_zeros(rsize,dom)); }
    VectorFunctionPatch<P,VARG> create_projection(VectorDomainType const& dom, Range indices) const {
        return VectorFunctionPatch<P,VARG>(this->_ptr->_create_projection(dom,indices)); }
    VectorFunctionPatch<P,VARG> create_identity(VectorDomainType const& dom) const {
        return VectorFunctionPatch<P,VARG>(this->_ptr->_create_identity(dom)); }

    ScalarFunctionPatch<P,SARG> create_zero(ScalarDomainType const& dom) const {
        return ScalarFunctionPatch<P,SARG>(this->_ptr->_create_zero(dom)); }
    VectorFunctionPatch<P,SARG> create_zeros(SizeType rsize, ScalarDomainType const& dom) const {
        return VectorFunctionPatch<P,SARG>(this->_ptr->_create_zeros(rsize,dom)); }
    ScalarFunctionPatch<P,SARG> create_constant(ScalarDomainType const& dom, Number<P> const& c) const {
        return ScalarFunctionPatch<P,SARG>(this->_ptr->_create_constant(dom,c)); }
    VectorFunctionPatch<P,SARG> create_constants(ScalarDomainType const& dom, Vector<Number<P>> const& c) const {
        return VectorFunctionPatch<P,SARG>(this->_ptr->_create_constants(dom,c)); }
    ScalarFunctionPatch<P,SARG> create_identity(ScalarDomainType const& dom) const {
        return ScalarFunctionPatch<P,SARG>(this->_ptr->_create_identity(dom)); }

    friend OutputStream& operator<<(OutputStream& os, FunctionPatchFactory<P> const& factory) { return factory._ptr->_write(os); }
};

template<class P> auto inline
FunctionPatchFactoryInterface<P>::create_identity(IntervalDomainType const& dom) const -> ScalarFunctionPatch<P,VARG> {
    return ScalarFunctionPatch<P,VARG>(this->create_coordinate(BoxDomainType(1u,dom),0u)); }

template<class FCTRY, class ARG> class FunctionPatchCreator {
    typedef typename FCTRY::Paradigm P;
  public:
    typedef FCTRY FactoryType;
    typedef typename DomainTraits<ARG>::BoundedDomainType DomainType;
    typedef P Paradigm;

    explicit FunctionPatchCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    decltype(auto) create(ScalarFunction<P,ARG> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create(VectorFunction<P,ARG> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create_zero() { return this->_factory.create_zero(this->_domain); }
    decltype(auto) create_zeros(SizeType n) { return this->_factory.create_zeros(n,this->_domain); }
    decltype(auto) create_constant(Number<P> const& c) const { return this->_factory.create_constant(this->_domain,c); }
    decltype(auto) create_identity() { return this->_factory.create_identity(this->_domain); }
  protected:
    FactoryType _factory;
    DomainType _domain;
};


//! \ingroup FunctionModule
//! \ingroup FunctionPatchSubModule
//! \brief Generic functions on bounded domains.
template<class P, class SIG> class FunctionPatch;

//! \ingroup FunctionPatchSubModule
//! \brief Generic scalar functions on bounded domains.
template<class P, class... ARGS> class FunctionPatch<P,RealScalar(ARGS...)>
//    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionPatch<P,D>, Number<P>>
    : public DispatchElementaryAlgebraOperations<ScalarFunctionPatch<P,ARGS...>, Number<P>>
{
    using RES=RealScalar; using SIG=RES(ARGS...);
  public:
    typedef Function<P,SIG> GenericType;
    typedef P Paradigm;
    typedef typename FunctionPatchInterface<P,SIG>::DomainType DomainType;
    typedef typename FunctionPatchInterface<P,SIG>::CodomainType CodomainType;
    typedef typename FunctionPatchInterface<P,SIG>::RangeType RangeType;
    typedef ExactNumber CoefficientType;
    typedef PositiveValidatedUpperNumber ErrorType;
    typedef Number<P> NumericType;
    typedef PositiveValidatedUpperNumber NormType;
    typedef typename ElementTraits<DomainType>::SizeType ArgumentSizeType;
    typedef typename ElementTraits<DomainType>::IndexType ArgumentIndexType;

    template<class Y> using Argument = typename ElementTraits<DomainType>::template Type<Y>;
    template<class Y> using Result = typename ElementTraits<CodomainType>::template Type<Y>;
  public:
    clone_on_copy_ptr< FunctionPatchInterface<P,SIG> > _ptr;
  public:
    FunctionPatch() : _ptr() { }
    explicit FunctionPatch(FunctionPatchInterface<P,SIG>* p) : _ptr(p) { }
    FunctionPatch(const SharedPointer<const FunctionPatchInterface<P,SIG>> p) : _ptr(p->_clone()) { }
    FunctionPatch(const FunctionPatch<P,SIG>& f) : _ptr(f._ptr) { }
    FunctionPatch& operator=(const FunctionPatch<P,SIG>& f) { this->_ptr=f._ptr; return *this; }
        FunctionPatch(const FunctionPatchInterface<P,SIG>& f) : _ptr(f._clone()) { }
    FunctionPatch(const Function<P,SIG>& f) : _ptr(dynamic_cast<FunctionPatchInterface<P,SIG>*>(f.raw_pointer()->_clone())) { }
    template<class PR, class PRE> inline explicit FunctionPatch(FunctionModel<P,SIG,PR,PRE> fm);
    operator Function<P,SIG>() const { return Function<P,SIG>(this->_ptr->_clone()); }
    operator FunctionPatchInterface<P,SIG>& () { return *_ptr; }
    operator const FunctionPatchInterface<P,SIG>& () const { return *_ptr; }
    const FunctionPatchInterface<P,SIG>* raw_pointer() const { return _ptr.operator->(); }
    FunctionPatchInterface<P,SIG>& reference() { return *_ptr; }
    const FunctionPatchInterface<P,SIG>& reference() const { return *_ptr; }

    ScalarFunctionPatch<P,ARGS...>& operator=(const Number<P>& c);
    ScalarFunctionPatch<P,ARGS...>& operator=(const ScalarFunction<P,ARGS...>& f);
    ScalarFunctionPatch<P,ARGS...>& operator=(const ScalarFunctionPatchInterface<P,ARGS...>& f);
//    ScalarFunctionPatch<P,D>& operator=(const ValidatedScalarMultivariateTaylorFunctionPatchDP& f);

    inline ArgumentSizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X operator() (const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class X> X evaluate(const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const; // { return this->_ptr->_range(); }

//    inline CoefficientType value() const { return this->_ptr->value(); }
//    inline CoefficientType gradient_value(ArgumentIndexType j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const; // { return this->_ptr->_error(); }
    inline Void clobber() const; // { return this->_ptr->_clobber(); }

//    inline ScalarFunctionPatch<P,D> apply(UnaryElementaryOperator op) const { return ScalarFunctionPatch<P,D>(this->_ptr->_apply(op)); }
    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionPatchCreator<FunctionPatchFactory<P>,ARGS...> factory(ScalarFunctionPatch<P,ARGS...> const& f) {
        FunctionPatchFactory<P> factory(f._ptr->_patch_factory());
        return FunctionPatchCreator<FunctionPatchFactory<P>,ARGS...>(f.domain(),factory); }
  public:
  public:
    friend Number<P> evaluate(const ScalarFunctionPatch<P,ARGS...>& f, const Vector<Number<P>>& x) {
        return f._ptr->_evaluate(x); }
    friend Number<P> unchecked_evaluate(const ScalarFunctionPatch<P,ARGS...>& f, const Vector<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(x); }

    friend ScalarFunctionPatch<P,ARGS...> partial_evaluate(const ScalarFunctionPatch<P,ARGS...>& f, SizeType j, const Number<P>& c) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_partial_evaluate(j,c)); }

    friend NormType norm(const ScalarFunctionPatch<P,ARGS...>& f) {
        return f._ptr->_generic_norm(); }
    friend ScalarFunctionPatch<P,ARGS...> derivative(const ScalarFunctionPatch<P,ARGS...>& f, SizeType j) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_derivative(j)); }
    friend ScalarFunctionPatch<P,ARGS...> antiderivative(const ScalarFunctionPatch<P,ARGS...>& f, SizeType j) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_antiderivative(j)); }
    friend ScalarFunctionPatch<P,ARGS...> antiderivative(const ScalarFunctionPatch<P,ARGS...>& f, SizeType j, Number<P> c) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_antiderivative(j,c)); }

    friend ScalarFunctionPatch<P,ARGS...> embed(const DomainType& d1, const ScalarFunctionPatch<P,ARGS...>& f, const DomainType& d2) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_embed(d1,d2)); }
    friend ScalarFunctionPatch<P,ARGS...> embed(const DomainType& d, const ScalarFunctionPatch<P,ARGS...>& f) {
        return embed(d,f,DomainType()); }
    friend ScalarFunctionPatch<P,ARGS...> embed(const ScalarFunctionPatch<P,ARGS...>& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend ScalarFunctionPatch<P,ARGS...> embed(const ScalarFunctionPatch<P,ARGS...>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend ScalarFunctionPatch<P,ARGS...> restrict(const ScalarFunctionPatch<P,ARGS...>& f, const DomainType& d) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_restriction(d)); }
    friend ScalarFunctionPatch<P,ARGS...> restriction(const ScalarFunctionPatch<P,ARGS...>& f, const DomainType& d) {
        return ScalarFunctionPatch<P,ARGS...>(f._ptr->_restriction(d)); }
    friend inline ScalarFunction<P,ARGS...> cast_unrestricted(ScalarFunctionPatch<P,ARGS...> const& f);

    friend VectorFunctionPatch<P,ARGS...> join(const ScalarFunctionPatch<P,ARGS...>& f1, const ScalarFunctionPatch<P,ARGS...>& f2) {
        return join(VectorFunctionPatch<P,ARGS...>(1,f1),f2); }
    friend VectorFunctionPatch<P,ARGS...> combine(const ScalarFunctionPatch<P,ARGS...>& f1, const ScalarFunctionPatch<P,ARGS...>& f2);
  public:
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionPatch<P,ARGS...>& f) {
        return os <<  f.operator ScalarMultivariateFunction<P>(); }
};

template<class P, class D> struct AlgebraOperations<ScalarFunctionPatch<P,D>> {

    typedef ScalarFunctionPatch<P,D> FM; typedef Number<P> X;
    typedef ScalarFunctionPatchInterface<P,D> FMI;

    static FM apply(BinaryElementaryOperator op, const FM& f1, const FM& f2) {
        return FM(&dynamic_cast<FMI&>(*f1._ptr->_apply(op,*f2._ptr))); }
    static FM apply(BinaryElementaryOperator op, const FM& f1, const X& c2) {
        return FM(&dynamic_cast<FMI&>(*f1._ptr->_apply(op,c2))); }
    static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
        return FM(&dynamic_cast<FMI&>(*f2._ptr->_rapply(op,c1))); }
    static FM apply(UnaryElementaryOperator op, const FM& f) {
        return FM(&dynamic_cast<FMI&>(*f._ptr->_apply(op))); }
    static FM apply(GradedElementaryOperator op, const FM& f, Int n) {
        return FM(&dynamic_cast<FMI&>(*f._ptr->_apply(op,n))); }
};


template<class P, class... ARGS> inline
ScalarFunctionPatch<P,ARGS...>& ScalarFunctionPatch<P,ARGS...>::operator=(const Number<P>& c) {
    return (*this)=factory(*this).create_constant(c); }
template<class P, class... ARGS> inline
ScalarFunctionPatch<P,ARGS...>& ScalarFunctionPatch<P,ARGS...>::operator=(const ScalarFunction<P,ARGS...>& f) {
    return (*this)=factory(*this).create(f); }
template<class P, class... ARGS> inline
ScalarFunctionPatch<P,ARGS...>& ScalarFunctionPatch<P,ARGS...>::operator=(const ScalarFunctionPatchInterface<P,ARGS...>& f) {
    return (*this)=ScalarFunctionPatch<P,ARGS...>(f._clone()); }



template<class V> struct Element;

template<class M> class ScaledFunctionPatch;
template<class M> class VectorScaledFunctionPatch;
//template<class M> struct Element<VectorScaledFunctionPatch<M>> { typedef ScaledFunctionPatch<M> Type; };

template<class P, class D> class VectorFunctionPatchElement
    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionPatch<P,D>, Number<P>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionPatch<P,D>>
{
    VectorFunctionPatch<P,D>* _p; SizeType _i;
  public:
    typedef typename ScalarFunctionPatch<P,D>::GenericType GenericType;
    operator const ScalarFunctionPatch<P,D> () const;
    VectorFunctionPatchElement(VectorFunctionPatch<P,D>* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionPatchElement<P,D>& operator=(const ScalarFunctionPatch<P,D>& sf) {
        _p->set(_i,sf); return *this; }
    VectorFunctionPatchElement<P,D>& operator=(const VectorFunctionPatchElement<P,D>& sf) {
        return this->operator=(static_cast<ScalarFunctionPatch<P,D>const>(sf)); }
    Void clobber() { ScalarFunctionPatch<P,D> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    const ScalarFunctionPatch<P,D> model() const { ScalarFunctionPatch<P,D> sf=_p->get(_i); return sf.model(); }
    const PositiveValidatedUpperNumber error() const { ScalarFunctionPatch<P,D> sf=_p->get(_i); return sf.error(); }
    friend inline OutputStream& operator<<(OutputStream& os, const VectorFunctionPatchElement<P,D>& function) {
        return os << static_cast< const ScalarFunctionPatch<P,D> >(function); }
};

//! \ingroup FunctionPatchSubModule
//! \brief Generic vector functions on bounded domains.
template<class P, class... ARGS> class FunctionPatch<P,RealVector(ARGS...)>
{
    using RES=RealVector; using SIG=RES(ARGS...);
  public:
    clone_on_copy_ptr< VectorFunctionPatchInterface<P,ARGS...> > _ptr;
  public:
    typedef P Paradigm;
    typedef VectorFunction<P,ARGS...> GenericType;
    typedef typename FunctionPatchInterface<P,SIG>::DomainType DomainType;
    typedef typename FunctionPatchInterface<P,SIG>::CodomainType CodomainType;
    typedef typename FunctionPatchInterface<P,SIG>::RangeType RangeType;
    typedef ExactNumber CoefficientType;
    typedef PositiveValidatedUpperNumber ErrorType;
    typedef Number<P> NumericType;
    typedef PositiveValidatedUpperNumber NormType;
    typedef typename ElementTraits<DomainType>::SizeType ArgumentSizeType;
    typedef typename ElementTraits<DomainType>::IndexType ArgumentIndexType;
    public:
    inline FunctionPatch() : _ptr() { }
    inline FunctionPatch(SharedPointer<const FunctionPatchInterface<P,SIG>> vfp)
        : _ptr(vfp->_clone()) { }
    inline FunctionPatch(SizeType n, const ScalarFunctionPatchInterface<P,ARGS...>& sf) {
        FunctionPatchFactory<P> factory(sf._patch_factory()); *this=factory.create_zeros(n,sf.domain());
        for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline FunctionPatch(Array<ScalarFunctionPatch<P,ARGS...>> const& asf)
        : FunctionPatch(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline FunctionPatch(List<ScalarFunctionPatch<P,ARGS...>> const& lsf)
        : FunctionPatch(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
    inline explicit FunctionPatch(FunctionPatchInterface<P,SIG>* p) : _ptr(p) { }
    template<class PR, class PRE> inline explicit FunctionPatch(FunctionModel<P,SIG,PR,PRE> fm);
    inline FunctionPatch(const FunctionPatchInterface<P,SIG>& f) : _ptr(f._clone()) { }
    inline FunctionPatch(const FunctionPatch<P,SIG>& f) : _ptr(f._ptr) { }
    inline FunctionPatch& operator=(const FunctionPatch<P,SIG>& f) { this->_ptr=f._ptr; return *this; }
    inline operator const FunctionPatchInterface<P,SIG>& () const { return *_ptr; }
    inline operator Function<P,SIG> () const { return Function<P,SIG>(*_ptr); }
    inline const FunctionPatchInterface<P,SIG>* raw_pointer() const { return _ptr.operator->(); }
    inline const FunctionPatchInterface<P,SIG>& reference() const { return *_ptr; }
    inline FunctionPatchInterface<P,SIG>& reference() { return *_ptr; }

    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    inline ScalarFunctionPatch<P,ARGS...> const get(SizeType i) const { return ScalarFunctionPatch<P,ARGS...>(this->_ptr->_get(i)); }
    inline Void set(SizeType i, ScalarFunctionPatch<P,ARGS...> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionPatch<P,ARGS...> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionPatchElement<P,ARGS...> operator[](SizeType i) { return VectorFunctionPatchElement<P,ARGS...>(this,i); }
    inline VectorFunctionPatch<P,ARGS...> operator[](Range rng) { VectorFunctionPatch<P,ARGS...> r=factory(*this).create_zeros(rng.size());
        for(SizeType i=0; i!=rng.size(); ++i) { r[i]=this->operator[](rng[i]); } return r; }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const; // { return this->_ptr->_range(); }
//        inline Vector<ErrorType> const errors() const { return this->_ptr->_errors(); }
//        inline ErrorType const error() const { return this->_ptr->_error(); }
    inline Void clobber(); // { this->_ptr->clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const;

    inline Void restrict(const DomainType& d) { *this=VectorFunctionPatch<P,ARGS...>(this->_ptr->_restriction(d)); }
//   public:
    friend FunctionPatchCreator<FunctionPatchFactory<P>,ARGS...> factory(VectorFunctionPatch<P,ARGS...> const& f) {
        FunctionPatchFactory<P> factory(f._ptr->_patch_factory());
        return FunctionPatchCreator<FunctionPatchFactory<P>,ARGS...>(f.domain(),factory); }
  public:
    friend inline ScalarFunctionPatch<P,ARGS...> compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionPatch<P,ARGS...>& g); // {
        // return ScalarFunctionPatch<P,ARGS...>(g._ptr->_compose(f)); }
    friend inline ScalarFunctionPatch<P,ARGS...> compose(const ScalarFunctionPatch<P,ARGS...>& f, const VectorFunctionPatch<P,ARGS...>& g) {
        return ScalarFunctionPatch<P,ARGS...>(g._ptr->_compose(f)); }
    friend inline VectorFunctionPatch<P,ARGS...> compose(const VectorMultivariateFunction<P>& f, const VectorFunctionPatch<P,ARGS...>& g); // {
        // return VectorFunctionPatch<P,ARGS...>(g._ptr->_compose(f)); }
    friend inline VectorFunctionPatch<P,ARGS...> compose(const VectorFunctionPatch<P,ARGS...>& f, const VectorFunctionPatch<P,ARGS...>& g) {
        return VectorFunctionPatch<P,ARGS...>(g._ptr->_compose(f)); }

    friend inline ScalarFunctionPatch<P,ARGS...> unchecked_compose(const ScalarFunctionPatch<P,ARGS...>& f, const VectorFunctionPatch<P,ARGS...>& g) {
        return ScalarFunctionPatch<P,ARGS...>(g._ptr->_unchecked_compose(f)); }
    friend inline VectorFunctionPatch<P,ARGS...> unchecked_compose(const VectorFunctionPatch<P,ARGS...>& f, const VectorFunctionPatch<P,ARGS...>& g) {
        return VectorFunctionPatch<P,ARGS...>(g._ptr->_unchecked_compose(f)); }

    friend inline VectorFunctionPatch<P,ARGS...> operator+(const VectorFunctionPatch<P,ARGS...>& f) {
        return VectorFunctionPatch<P,ARGS...>(f._ptr->_clone()); }
    friend inline VectorFunctionPatch<P,ARGS...> operator-(const VectorFunctionPatch<P,ARGS...>& f) {
        VectorFunctionPatch<P,ARGS...> r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionPatch<P,ARGS...> operator+(const VectorFunctionPatch<P,ARGS...>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        VectorFunctionPatch<P,ARGS...> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionPatch<P,ARGS...> operator-(const VectorFunctionPatch<P,ARGS...>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        VectorFunctionPatch<P,ARGS...> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionPatch<P,ARGS...> operator+(const VectorFunctionPatch<P,ARGS...>& f1, const Vector<Number<P>>& c2) {
        VectorFunctionPatch<P,ARGS...> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionPatch<P,ARGS...> operator-(const VectorFunctionPatch<P,ARGS...>& f1, const Vector<Number<P>>& c2) {
        VectorFunctionPatch<P,ARGS...> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionPatch<P,ARGS...> operator+(const Vector<Number<P>>& c1, const VectorFunctionPatch<P,ARGS...>& f2);
    friend inline VectorFunctionPatch<P,ARGS...> operator-(const Vector<Number<P>>& c1, const VectorFunctionPatch<P,ARGS...>& f2);
    friend inline VectorFunctionPatch<P,ARGS...> operator*(const VectorFunctionPatch<P,ARGS...>& f1, const Number<P>& c2) {
        VectorFunctionPatch<P,ARGS...> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionPatch<P,ARGS...> operator*(const Number<P>& c1, const VectorFunctionPatch<P,ARGS...>& f2) {
        VectorFunctionPatch<P,ARGS...> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

    friend inline VectorFunctionPatch<P,ARGS...> operator+(const VectorFunctionPatch<P,ARGS...>& f1, const VectorMultivariateFunction<P>& f2) {
        return f1+factory(f1).create(f2); }
    friend inline VectorFunctionPatch<P,ARGS...> operator-(const VectorFunctionPatch<P,ARGS...>& f1, const VectorMultivariateFunction<P>& f2) {
        return f1-factory(f1).create(f2); }
    friend inline VectorFunctionPatch<P,ARGS...> operator*(const VectorFunctionPatch<P,ARGS...>& f1, const ScalarMultivariateFunction<P>& f2) {
        return f1*factory(f1).create(f2); }
    friend inline VectorFunctionPatch<P,ARGS...> operator/(const VectorFunctionPatch<P,ARGS...>& f1, const ScalarMultivariateFunction<P>& f2) {
        return f1/factory(f1).create(f2); }
    friend inline VectorFunctionPatch<P,ARGS...> operator+(const VectorMultivariateFunction<P>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return factory(f2).create(f1)+f2; }
    friend inline VectorFunctionPatch<P,ARGS...> operator-(const VectorMultivariateFunction<P>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return factory(f2).create(f1)-f2; }
    friend inline VectorFunctionPatch<P,ARGS...> operator*(const ScalarMultivariateFunction<P>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return factory(f2).create(f1)*f2; }
    friend inline VectorFunctionPatch<P,ARGS...> operator/(const ScalarMultivariateFunction<P>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return factory(f2).create(f1)/f2; }


  public:
    friend NormType norm(const VectorFunctionPatch<P,ARGS...>& f) {
        return f._ptr->_generic_norm(); }
    friend VectorFunctionPatch<P,ARGS...> embed(const DomainType& d1, const VectorFunctionPatch<P,ARGS...>& f, const DomainType& d2) {
        return VectorFunctionPatch<P,ARGS...>(f._ptr->_embed(d1,d2)); }
    friend VectorFunctionPatch<P,ARGS...> embed(const DomainType& d, const VectorFunctionPatch<P,ARGS...>& f) {
        return embed(d,f,DomainType()); }
    friend VectorFunctionPatch<P,ARGS...> embed(const VectorFunctionPatch<P,ARGS...>& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend VectorFunctionPatch<P,ARGS...> embed(const VectorFunctionPatch<P,ARGS...>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend VectorFunctionPatch<P,ARGS...> restrict(const VectorFunctionPatch<P,ARGS...>& f, const DomainType& d) {
        VectorFunctionPatchInterface<P,ARGS...>* rptr=f._ptr->_clone(); rptr->restrict(d); return VectorFunctionPatch<P,ARGS...>(rptr); }
    friend VectorFunctionPatch<P,ARGS...> restriction(const VectorFunctionPatch<P,ARGS...>& f, const DomainType& d) {
        VectorFunctionPatchInterface<P,ARGS...>* rptr=f._ptr->_clone(); rptr->restrict(d); return VectorFunctionPatch<P,ARGS...>(rptr); }
    friend inline VectorFunction<P,ARGS...> cast_unrestricted(VectorFunctionPatch<P,ARGS...> const& f);

    friend Vector<Number<P>> evaluate(const VectorFunctionPatch<P,ARGS...>& f, const Vector<Number<P>>& x) {
        return f._ptr->_evaluate(x); }

    friend Vector<Number<P>> unchecked_evaluate(const VectorFunctionPatch<P,ARGS...>& f, const Vector<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(x); }

    friend VectorFunctionPatch<P,ARGS...> partial_evaluate(const VectorFunctionPatch<P,ARGS...>& f, SizeType j, const Number<P>& c) {
        return VectorFunctionPatch<P,ARGS...>(f._ptr->_partial_evaluate(j,c)); }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionPatch<P,ARGS...>& f) {
        return os <<  f.operator VectorMultivariateFunction<P>(); }

    friend ScalarFunctionPatch<P,ARGS...> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionPatch<P,ARGS...>& g) {
        ScalarFunctionPatchInterface<P,ARGS...> const* fptr = dynamic_cast<ScalarFunctionPatchInterface<P,ARGS...> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(ScalarFunctionPatch<P,ARGS...>(*fptr),g); } else { return compose(f,g); } }
    friend VectorFunctionPatch<P,ARGS...> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionPatch<P,ARGS...>& g) {
        VectorFunctionPatchInterface<P,ARGS...> const* fptr = dynamic_cast<VectorFunctionPatchInterface<P,ARGS...> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(VectorFunctionPatch<P,ARGS...>(*fptr),g); } else { return compose(f,g); } }

    friend VectorFunctionPatch<P,ARGS...> antiderivative(const VectorFunctionPatch<P,ARGS...>& f, SizeType j) {
        VectorFunctionPatch<P,ARGS...> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
        return r;
    }

    friend VectorFunctionPatch<P,ARGS...> antiderivative(const VectorFunctionPatch<P,ARGS...>& f, SizeType j, const Number<P>& c) {
        VectorFunctionPatch<P,ARGS...> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
        return r;
    }

    //friend VectorFunctionPatch<P,ARGS...> join(const ScalarFunctionPatch<P,ARGS...>& f1, const ScalarFunctionPatch<P,ARGS...>& f2);
    friend VectorFunctionPatch<P,ARGS...> join(const ScalarFunctionPatch<P,ARGS...>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return join(VectorFunctionPatch<P,ARGS...>(1u,f1),f2); }
    friend VectorFunctionPatch<P,ARGS...> join(const VectorFunctionPatch<P,ARGS...>& f1, const ScalarFunctionPatch<P,ARGS...>& f2) {
        VectorFunctionPatch<P,ARGS...> r=VectorFunctionPatch<P,ARGS...>(f1._ptr->_clone()); r._ptr->_adjoin(f2); return r; }
    friend VectorFunctionPatch<P,ARGS...> join(const VectorFunctionPatch<P,ARGS...>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return VectorFunctionPatch<P,ARGS...>(f1._ptr->_join(f2)); }

    friend VectorFunctionPatch<P,ARGS...> combine(const ScalarFunctionPatch<P,ARGS...>& f1, const ScalarFunctionPatch<P,ARGS...>& f2) {
        return VectorFunctionPatch<P,ARGS...>(1,f1)._ptr->_combine(VectorFunctionPatch<P,ARGS...>(1,f2)); };
    friend VectorFunctionPatch<P,ARGS...> combine(const ScalarFunctionPatch<P,ARGS...>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return VectorFunctionPatch<P,ARGS...>(1,f1)._ptr->_combine(f2); };
    friend VectorFunctionPatch<P,ARGS...> combine(const VectorFunctionPatch<P,ARGS...>& f1, const ScalarFunctionPatch<P,ARGS...>& f2) {
        return VectorFunctionPatch<P,ARGS...>(f1._ptr->_combine(VectorFunctionPatch<P,ARGS...>(1,f2))); };
    friend VectorFunctionPatch<P,ARGS...> combine(const VectorFunctionPatch<P,ARGS...>& f1, const VectorFunctionPatch<P,ARGS...>& f2) {
        return VectorFunctionPatch<P,ARGS...>(f1._ptr->_combine(f2)); }
};

template<class P> inline Number<P> unchecked_evaluate(const ScalarMultivariateFunction<P>& f, const Vector<Number<P>>& x) {
    ScalarMultivariateFunctionPatchInterface<P> const* fptr = dynamic_cast<ScalarMultivariateFunctionPatchInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarMultivariateFunctionPatch<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P> inline Vector<Number<P>> unchecked_evaluate(const VectorMultivariateFunction<P>& f, const Vector<Number<P>>& x) {
    VectorMultivariateFunctionPatchInterface<P> const* fptr = dynamic_cast<VectorMultivariateFunctionPatchInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorMultivariateFunctionPatch<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P, class D> inline ScalarFunctionPatch<P,D> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionPatch<P,D>& g) {
    ScalarFunctionPatchInterface<P,D> const* fptr = dynamic_cast<ScalarFunctionPatchInterface<P,D> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionPatch<P,D>(*fptr),g); } else { return compose(f,g); } }
template<class P, class D> inline VectorFunctionPatch<P,D> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionPatch<P,D>& g) {
    VectorFunctionPatchInterface<P,D> const* fptr = dynamic_cast<VectorFunctionPatchInterface<P,D> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionPatch<P,D>(*fptr),g); } else { return compose(f,g); } }

// FIXME: Should be unneeded
template<class D> ScalarFunctionPatch<ValidatedTag,D> unchecked_compose(const ScalarMultivariateFunctionPatch<ValidatedTag>& f, const VectorFunctionPatch<ValidatedTag,D>& g) {
    return VectorFunctionPatch<ValidatedTag,D>(g._ptr->_unchecked_compose(f)); }
template<class D> VectorFunctionPatch<ValidatedTag,D> unchecked_compose(const VectorMultivariateFunctionPatch<ValidatedTag>& f, const VectorFunctionPatch<ValidatedTag,D>& g) {
    return VectorFunctionPatch<ValidatedTag,D>(g._ptr->_unchecked_compose(f)); }

// FIXME: Should be unneeded
template<class D> ScalarFunctionPatch<ValidatedTag,D> restrict(const ScalarFunctionPatch<ValidatedTag,D>& f, const D& dom) {
    return ScalarFunctionPatch<ValidatedTag,D>(f._ptr->_restriction(dom)); }
template<class D> VectorFunctionPatch<ValidatedTag,D> restrict(const VectorFunctionPatch<ValidatedTag,D>& f, const D& dom) {
    return VectorFunctionPatch<ValidatedTag,D>(f._ptr->_restriction(dom)); }


template<class P, class D> VectorFunctionPatchElement<P,D>::operator const ScalarFunctionPatch<P,D> () const {
    return ScalarFunctionPatch<P,D>(_p->get(_i)); }


} // namespace Ariadne

#endif
