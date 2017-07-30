/***************************************************************************
 *            symbolic_function.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file symbolic_function.hpp
 *  \brief Symbolic functions
 */

#ifndef ARIADNE_SYMBOLIC_FUNCTION_HPP
#define ARIADNE_SYMBOLIC_FUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function_interface.hpp"

#include "utility/macros.hpp"
#include "utility/pointer.hpp"
#include "utility/container.hpp"
#include "utility/metaprogramming.hpp"

#include "numeric/numeric.hpp"
#include "numeric/operators.tpl.hpp"
#include "algebra/vector.hpp"

#include "function/function_mixin.hpp"
#include "function/projection.hpp"
#include "function/formula.hpp"

namespace Ariadne {

//------------------------ Formula functions  -----------------------------------//

//! A function defined by a formula
template<class Y>
struct ScalarFormulaFunction
    : ScalarFunctionMixin<ScalarFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    SizeType _argument_size;
    Formula<Y> _formula;

    ScalarFormulaFunction(SizeType as, const Formula<Y>& f) : _argument_size(as), _formula(f) { }
    operator Formula<Y>() const { return _formula; }

    virtual SizeType argument_size() const final { return _argument_size; }
    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const final { return new ScalarFormulaFunction<Y>(_argument_size,Ariadne::derivative(_formula,j)); }
    virtual OutputStream& write(OutputStream& os) const final { return os << this->_formula; }
    virtual OutputStream& repr(OutputStream& os) const final { return os << "FormulaFunction("<<this->_argument_size<<","<<this->_formula<<")"; }
    template<class X> Void _compute(X& r, const Vector<X>& x) const { r=Ariadne::cached_evaluate(_formula,x); }
};

typedef ScalarFormulaFunction<EffectiveNumber> EffectiveScalarFormulaFunction;

//! A vector function defined by formulae
template<class Y>
struct VectorFormulaFunction
    : VectorFunctionMixin<VectorFormulaFunction<Y>,InformationTag<Y>>
{
    SizeType _argument_size;
    Vector< Formula<Y> > _formulae;

    VectorFormulaFunction(SizeType as, const List< Formula<Y> >& f) : _argument_size(as), _formulae(f) { }
    VectorFormulaFunction(SizeType as, const Vector< Formula<Y> >& f) : _argument_size(as), _formulae(f) { }

    virtual SizeType result_size() const { return this->_formulae.size(); }
    virtual SizeType argument_size() const { return this->_argument_size; }
    virtual ScalarFormulaFunction<Y>* _get(SizeType i) const { return new ScalarFormulaFunction<Y>(this->_argument_size,this->_formulae[i]); }
    virtual VectorFormulaFunction<Y>* _derivative(SizeType k) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& write(OutputStream& os) const { return os << this->_formulae; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "VectorFormulaFunction("<<this->result_size()<<","<<this->argument_size()<<","<<this->_formulae<<")"; }
    template<class X> Void _compute(Vector<X>& r, const Vector<X>& x) const { r=Ariadne::cached_evaluate(this->_formulae,x); }
};



//------------------------ Arithmetic scalar functions  -----------------------------------//


//! A constant function f(x)=c
template<class Y>
struct ConstantFunction
    : ScalarFunctionMixin<ConstantFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
  public:
    BoxDomainType _domain;
    Y _value;

    //ConstantFunction(SizeType as, const Y& c) : _argument_size(as), _value(c) { }
    ConstantFunction(BoxDomainType dom, const Y& c) : _domain(dom), _value(c) { }
    operator Y() const { return _value; }

    virtual const BoxDomainType domain() const { return _domain; }
    virtual SizeType argument_size() const { return _domain.dimension(); }
    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const { return new ConstantFunction<Y>(_domain,Y(0)); }
    virtual OutputStream& write(OutputStream& os) const { return os << this->_value; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "CF[R"<<this->argument_size()<<"]("<<_value<<")"; }
    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=make_constant(_value,x.zero_element()); }
};


//! A coordinate function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=x_i\f$.
template<class P>
struct CoordinateFunction
    : ScalarFunctionMixin<CoordinateFunction<P>,P>
{
    typedef Number<P> Y;

    BoxDomainType _domain;
    SizeType _index;

    //CoordinateFunction(SizeType as, SizeType i) : _argument_size(as), _index(i) { }
    CoordinateFunction(BoxDomainType dom, SizeType i) : _domain(dom), _index(i) { }
    SizeType index() const { return _index; }

    virtual const BoxDomainType domain() const { return _domain; }
    virtual SizeType argument_size() const { return _domain.dimension(); }
    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const {
        if(j==_index) { return new ConstantFunction<Y>(_domain,Y(1)); }
        else { return new ConstantFunction<Y>(_domain,Y(0)); } }
    virtual OutputStream& write(OutputStream& os) const { return os << "x"<<this->_index; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "IF[R"<<this->argument_size()<<"](x"<<this->_index<<")"; }
    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=x[_index]; }
};


//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
template<class P>
struct UnaryFunction
    : ScalarFunctionMixin< UnaryFunction<P>, P >
{
    typedef Number<P> Y;
  public:
    UnaryFunction(const OperatorCode& op, const ScalarFunction<P>& arg)
        : _op(op), _arg(arg) { }
    virtual UnaryFunction<P>* clone() const { return new UnaryFunction<P>(*this); }
    virtual const BoxDomainType domain() const {
        return this->_arg.domain(); }
    virtual SizeType argument_size() const {
        return this->_arg.argument_size(); }

    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const {
        return static_cast<const ScalarFunctionInterface<P>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<P> derivative(SizeType j) const {
        switch(_op) {
            case OperatorCode::POS: return _arg.derivative(j);
            case OperatorCode::NEG: return -_arg.derivative(j);
            case OperatorCode::REC: return -_arg.derivative(j)/sqr(_arg);
            case OperatorCode::SQR: return 2*_arg.derivative(j)*_arg;
            case OperatorCode::SQRT: return _arg.derivative(j)/(2*sqrt(_arg));
            case OperatorCode::EXP: return _arg*_arg.derivative(j);
            case OperatorCode::LOG: return _arg.derivative(j)/_arg;
            case OperatorCode::SIN: return _arg.derivative(j)*cos(_arg);
            case OperatorCode::COS: return -_arg.derivative(j)*sin(_arg);
            default: ARIADNE_FAIL_MSG("Unknown unary function "<<this->_op);
        }
    }

    virtual OutputStream& repr(OutputStream& os) const {
        return os << "UF[R" << this->argument_size() << "](" << *this << ")"; }
    virtual OutputStream& write(OutputStream& os) const {
        return os << _op << '(' << _arg << ')'; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=Ariadne::compute(_op,_arg.evaluate(x)); }

    OperatorCode _op;
    ScalarFunction<P> _arg;
};


template<class P> ScalarFunction<P> sqr(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(OperatorCode::SQR,f)); }
template<class P> ScalarFunction<P> sqrt(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(OperatorCode::SQRT,f)); }
template<class P> ScalarFunction<P> sin(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(OperatorCode::SIN,f)); }
template<class P> ScalarFunction<P> cos(const ScalarFunction<P>& f) {
    return ScalarFunction<P>(new UnaryFunction<P>(OperatorCode::COS,f)); }

template<class P>
struct BinaryFunction
    : ScalarFunctionMixin< BinaryFunction<P>, P >
{
    typedef Number<P> Y;
  public:
    BinaryFunction(OperatorCode op, const ScalarFunction<P>& arg1, const ScalarFunction<P>& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) { ARIADNE_ASSERT_MSG(arg1.argument_size()==arg2.argument_size(),"op='"<<op<<"', arg1="<<arg1<<", arg2="<<arg2); }
    virtual BinaryFunction<P>* clone() const { return new BinaryFunction<P>(*this); }
    virtual const BoxDomainType domain() const {
        return intersection(this->_arg1.domain(),this->_arg2.domain()); }
    virtual SizeType argument_size() const {
        return this->_arg1.argument_size(); }

    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const {
        return static_cast<const ScalarFunctionInterface<P>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<P> derivative(SizeType j) const {
        switch(_op) {
            case OperatorCode::ADD:
                return _arg1.derivative(j)+_arg2.derivative(j);
            case OperatorCode::SUB:
                return _arg1.derivative(j)-_arg2.derivative(j);
            case OperatorCode::MUL:
                return _arg1.derivative(j)*_arg2+_arg1*_arg2.derivative(j);
            case OperatorCode::DIV:
                if(dynamic_cast<const ConstantFunction<Y>*>(_arg2.raw_pointer())) {
                    return _arg1.derivative(j)/_arg2;
                } else {
                    return _arg1.derivative(j)/_arg2-_arg2.derivative(j)*_arg1/sqr(_arg2);
                }
            default: ARIADNE_FAIL_MSG("Unknown binary function "<<this->_op);
        }
    }

    virtual OutputStream& repr(OutputStream& os) const {
        return os << "BF[R" << this->argument_size() << "](" << *this << ")"; }
    virtual OutputStream& write(OutputStream& os) const {
        if(_op==OperatorCode::ADD || _op==OperatorCode::SUB) { return os << '(' << _arg1 << symbol(_op) << _arg2 << ')'; }
        else { return os << _arg1 << symbol(_op) << _arg2; } }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=Ariadne::compute(_op,_arg1.evaluate(x),_arg2.evaluate(x)); }

    OperatorCode _op;
    ScalarFunction<P> _arg1;
    ScalarFunction<P> _arg2;
};


// \brief The power function \f$(x,n)\mapsto x^n\f$.
template<class P>
class GradedFunction
    : public ScalarFunctionMixin< GradedFunction<P>, P >
{
    typedef Number<P> Y;
  public:
    GradedFunction(OperatorCode op, const ScalarFunction<P>& arg1, const Int& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {  }
    virtual GradedFunction<P>* clone() const { return new GradedFunction<P>(*this); }
    virtual const BoxDomainType domain() const {
        return this->_arg1.domain(); }
    virtual SizeType argument_size() const {
        return this->_arg1.argument_size(); }

    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const {
        ARIADNE_NOT_IMPLEMENTED;
    }

    virtual ScalarFunction<P> derivative(SizeType j) const {
        assert(_op==OperatorCode::POW);
        if(_arg2==0) { return ScalarFunction<P>::constant(this->argument_size(),Y(0)); }
        if(_arg2==1) { return _arg1.derivative(j); }
        if(_arg2==2) { return 2*_arg1.derivative(j)*_arg1; }
        return _arg2*_arg1.derivative(j)*pow(_arg1,_arg2-1);
    }

    virtual OutputStream& repr(OutputStream& os) const {
        return os << "GF["<<this->argument_size()<<"]("<< *this <<")"; }
    virtual OutputStream& write(OutputStream& os) const {
        return os << _op << "(" << _arg1 << "," << _arg2 << ")"; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=compute(_op,_arg1.evaluate(x),_arg2); }

    OperatorCode _op;
    ScalarFunction<P> _arg1;
    Int _arg2;
};

typedef ConstantFunction<Real> RealConstantFunction;
typedef ConstantFunction<EffectiveNumericType> EffectiveConstantFunction;
typedef CoordinateFunction<EffectiveTag> EffectiveCoordinateFunction;
typedef UnaryFunction<EffectiveTag> EffectiveUnaryFunction;
typedef BinaryFunction<EffectiveTag> EffectiveBinaryFunction;
typedef GradedFunction<EffectiveTag> EffectiveGradedFunction;


//------------------------ Vector of Scalar functions  -----------------------------------//

template<class P, class D=BoxDomainType> class NonResizableScalarFunction : public ScalarFunction<P,D> {
  public:
    NonResizableScalarFunction<P,D>& operator=(const ScalarFunction<P,D>& f) {
        ARIADNE_ASSERT_MSG(this->domain()==f.domain(), "this->domain()="<<this->domain()<<", f.domain()="<<f.domain()<<"\n\n*this="<<*this<<"\nf="<<f<<"\n\n");
        this->ScalarFunction<P,D>::operator=(f);
        return *this;
    }
};

template<class P, class D=BoxDomainType>
struct VectorOfScalarFunction
    : VectorFunctionMixin<VectorOfScalarFunction<P,D>,P,D>
    , public virtual VectorOfFunctionInterface<P,D>
{
    typedef D DomainType;
    VectorOfScalarFunction(SizeType rs, SizeType as)
        : VectorOfScalarFunction(rs, ScalarFunction<P,D>(as)) { }
    VectorOfScalarFunction(SizeType rs, DomainType dom)
        : VectorOfScalarFunction(rs, ScalarFunction<P,D>(dom)) { }
    VectorOfScalarFunction(SizeType rs, const ScalarFunction<P,D>& f)
        : _dom(f.domain()), _vec(rs,f) { }
    VectorOfScalarFunction(const Vector<ScalarFunction<P,D>>& vsf)
        : _dom(vsf.zero_element().domain()), _vec(vsf) { }

    Void set(SizeType i, ScalarFunction<P,D> f) {
        if(this->argument_size()==0u) { this->_dom=f.domain(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    ScalarFunction<P,D> get(SizeType i) const {
        return this->_vec[i]; }

    virtual SizeType result_size() const final {
        return _vec.size(); }
    virtual SizeType argument_size() const final {
        return _dom.dimension(); }
    virtual DomainType const domain() const final {
        return _dom; }

    virtual ScalarFunctionInterface<P,D>* _get(SizeType i) const final {
        return this->_vec[i].raw_pointer()->_clone(); }
    virtual Void _set(SizeType i, const ScalarFunctionInterface<P,D>* sf) final {
        this->_vec[i]=ScalarFunction<P,D>(sf->_clone()); }
    virtual VectorFunctionInterface<P,D>* _derivative(SizeType i) const {
        ARIADNE_NOT_IMPLEMENTED; }

    const ScalarFunction<P,D> operator[](SizeType i) const {
        return this->_vec[i]; }

    NonResizableScalarFunction<P,D>& operator[](SizeType i) {
        return static_cast<NonResizableScalarFunction<P,D>&>(this->_vec[i]); }

    virtual OutputStream& write(OutputStream& os) const {
        os << "[";
        for(SizeType i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->write(os); }
        return os << "]"; }

    virtual OutputStream& repr(OutputStream& os) const {
        //os << "VoSF[R" << this->argument_size() << "->R" << this->result_size() << "]";
        os << "[";
        for(SizeType i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->repr(os); }
        return os << "]"; }

    template<class X> inline Void _compute(Vector<X>& r, const ElementType<D,X>& x) const {
        r=Vector<X>(this->_vec.size(),zero_element(x));
        for(SizeType i=0; i!=r.size(); ++i) {
            r[i]=_vec[i].evaluate(x); } }

    DomainType _dom;
    Vector<ScalarFunction<P,D>> _vec;

};


template<class P, class D=BoxDomainType>
struct FunctionElement
    : ScalarFunctionMixin<FunctionElement<P,D>,P>
{
    typedef D DomainType;

    FunctionElement(const VectorFunction<P,D>& f, SizeType i)
        : _f(f), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }

    virtual SizeType argument_size() const { return _f.argument_size(); }
    virtual DomainType domain() const { return _f.domain(); }
    virtual OutputStream& write(OutputStream& os) const { return os<<_f<<"["<<_i<<"]"; }
    virtual ScalarFunctionInterface<P,D>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=this->_f.evaluate(x)[_i]; }

    VectorFunction<P,D> _f;
    SizeType _i;
};

//------------------------ Results of functional operations  -----------------------------------//

template<class P>
struct ScalarEmbeddedFunction
    : ScalarFunctionMixin<ScalarEmbeddedFunction<P>,P>
{
    ScalarEmbeddedFunction(SizeType as1, const ScalarFunction<P>& f2, SizeType as3)
        : _as1(as1), _f2(f2), _as3(as3) { }
    virtual SizeType argument_size() const { return _as1+_f2.argument_size()+_as3; }
    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& write(OutputStream& os) const { return os << "ScalarEmbeddedFunction( as1="<<_as1<<", f2="<<_f2<<", as3="<<_as3<<" )"; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        Vector<X> px=project(x,Range(_as1,_as1+_f2.argument_size())); r=_f2.evaluate(px); }

    SizeType _as1;
    ScalarFunction<P> _f2;
    SizeType _as3;
};


template<class P>
struct VectorEmbeddedFunction
    : VectorFunctionMixin<VectorEmbeddedFunction<P>,P>
{
    VectorEmbeddedFunction(SizeType as1, const VectorFunction<P>& f2, SizeType as3)
        : _as1(as1), _f2(f2), _as3(as3) { }
    virtual SizeType result_size() const { return _f2.result_size(); }
    virtual SizeType argument_size() const { return _as1+_f2.argument_size()+_as3; }
    virtual ScalarFunctionInterface<P>* _get(SizeType i) const { return new ScalarEmbeddedFunction<P>(_as1,_f2.get(i),_as3); }
    virtual VectorFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& write(OutputStream& os) const { return os << "VectorEmbeddedFunction( as1="<<_as1<<", f2="<<_f2<<", as3="<<_as3<<" )"; }

    template<class X> inline Void _compute(Vector<X>& r, const Vector<X>& x) const {
        Vector<X> px=project(x,Range(_as1,_as1+_f2.argument_size())); r=_f2.evaluate(px); }

    SizeType _as1;
    VectorFunction<P> _f2;
    SizeType _as3;
};


template<class P>
struct ScalarComposedFunction
    : ScalarFunctionMixin<ScalarComposedFunction<P>,P>
{
    ScalarComposedFunction(const ScalarFunction<P>& f, const VectorFunction<P>& g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual SizeType argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& write(OutputStream& os) const { return os << "ScalarComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    ScalarFunction<P> _f;
    VectorFunction<P> _g;
};


template<class P>
struct VectorComposedFunction
    :  VectorFunctionMixin<VectorComposedFunction<P>,P>
{
    VectorComposedFunction(VectorFunction<P> f, VectorFunction<P> g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual SizeType result_size() const { return _f.result_size(); }
    virtual SizeType argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<P>* _get(SizeType i) const { return new ScalarComposedFunction<P>(_f[i],_g); }
    virtual VectorFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& write(OutputStream& os) const { return os << "ComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class X> inline Void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    VectorFunction<P> _f;
    VectorFunction<P> _g;
};


template<class P>
struct JoinedFunction
    : VectorFunctionMixin<JoinedFunction<P>,P>
{
    typedef BoxDomainType D;
    JoinedFunction(VectorFunction<P,D> f1, VectorFunction<P,D> f2)
        : _f1(f1), _f2(f2) { ARIADNE_ASSERT(f1.argument_size()==f2.argument_size()); }
    virtual SizeType result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual SizeType argument_size() const { return _f1.argument_size(); }
    virtual OutputStream& write(OutputStream& os) const { return os << "JoinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }
    virtual ScalarFunctionInterface<P,D>* _get(SizeType i) const {
        return (i<_f1.result_size()) ? dynamic_cast<VectorOfFunctionInterface<P,D>const*>(_f1.raw_pointer())->_get(i)
                                     : dynamic_cast<VectorOfFunctionInterface<P,D>const*>(_f2.raw_pointer())->_get(i-_f1.result_size()); }
    virtual VectorFunctionInterface<P,D>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    template<class X> inline Void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=join(_f1.evaluate(x),_f2.evaluate(x)); }

    VectorFunction<P,D> _f1;
    VectorFunction<P,D> _f2;
};


template<class P>
class CombinedFunction
    : VectorFunctionMixin<CombinedFunction<P>,P>
{
    CombinedFunction(VectorFunction<P> f1, VectorFunction<P> f2)
        : _f1(f1), _f2(f2) { }
    virtual SizeType result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual SizeType argument_size() const { return _f1.argument_size()+_f2.argument_size(); }
    virtual OutputStream& write(OutputStream& os) const { return os << "CombinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }

    virtual VectorFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> inline Void _compute(Vector<X>& r, const Vector<X>& x) const {
        return r=combine(_f1.evaluate(project(x,range(0,_f1.argument_size()))),
                         _f2.evaluate(project(x,range(_f1.argument_size(),this->argument_size())))); }

    VectorFunction<P> _f1;
    VectorFunction<P> _f2;
};


template<class P>
class ProjectedFunction
    : VectorFunctionMixin<ProjectedFunction<P>,P>
{
    ProjectedFunction(VectorFunction<P> f, Projection prj)
        : _f(f), _prj(prj) { ARIADNE_PRECONDITION(f.result_size()==prj.argument_size()); }
    virtual SizeType result_size() const { return _prj.result_size(); }
    virtual SizeType argument_size() const { return _f.argument_size(); }
    virtual OutputStream& write(OutputStream& os) const { return os << "ProjectedFunction( f="<<_f<<", prj="<<_prj<<" )"; }

    virtual VectorFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> inline Void _compute(Vector<X>& r, const Vector<X>& x) const {
        return r=_prj(_f(x)); }
    VectorFunction<P> _f;
    Projection _prj;
};


// A Lie deriviative \f$\nabla g\cdot f\f$.
template<class P>
struct LieDerivativeFunction
    : ScalarFunctionMixin<LieDerivativeFunction<P>,P>
{
    //! \brief Construct the identity function in dimension \a n.
    LieDerivativeFunction(const ScalarFunction<P>& g, const VectorFunction<P>& f) {
        ARIADNE_ASSERT(g.argument_size()==f.argument_size());
        ARIADNE_ASSERT(f.result_size()==f.argument_size());
        _g=g; for(SizeType j=0; j!=g.argument_size(); ++j) { _dg[j]=g.derivative(j); } _f=f; }
    SizeType argument_size() const { return _g.argument_size(); }
    virtual ScalarFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& write(OutputStream& os) const { return os << "LieDerivative( g="<<_g<<", f="<<_f<<" )"; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        //const Vector<R> fx=_f.evaluate(x); r=0; for(SizeType i=0; i!=_dg.size(); ++i) { r+=fx[i]+_dg[i].evaluate(x); } }
        Vector<X> fx=_f.evaluate(x);
        r=0;
        for(SizeType i=0; i!=_dg.size(); ++i) {
            r+=fx[i]+_dg[i].evaluate(x);
        }
    }

    ScalarFunction<P> _g;
    List< ScalarFunction<P> > _dg;
    VectorFunction<P> _f;
};




} // namespace Ariadne

#endif