/***************************************************************************
 *            function/symbolic_function.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file function/symbolic_function.hpp
 *  \brief Symbolic functions
 */

#ifndef ARIADNE_SYMBOLIC_FUNCTION_HPP
#define ARIADNE_SYMBOLIC_FUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../function/function_interface.hpp"

#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"
#include "../utility/metaprogramming.hpp"

#include "../numeric/numeric.hpp"
#include "../numeric/operators.tpl.hpp"
#include "../algebra/vector.hpp"

#include "../function/function_mixin.hpp"
#include "../function/projection.hpp"
#include "../function/formula.hpp"

namespace Ariadne {

//------------------------ Formula functions  -----------------------------------//

//! A function defined by a formula
template<class Y>
struct ScalarUnivariateFormulaFunction
    : ScalarUnivariateFunctionMixin<ScalarUnivariateFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    Formula<Y> _formula;

    ScalarUnivariateFormulaFunction(const Formula<Y>& f) : _formula(f) { }
    operator Formula<Y>() const { return _formula; }

    virtual SizeOne argument_size() const final { return SizeOne(); }
    virtual SizeOne result_size() const final { return SizeOne(); }
    virtual typename ScalarUnivariateFunction<P>::Interface* _derivative(SizeOne j) const final {
        return new ScalarUnivariateFormulaFunction<Y>(Ariadne::derivative(_formula,0)); }
    virtual OutputStream& _write(OutputStream& os) const final { return os << this->_formula; }
    virtual OutputStream& repr(OutputStream& os) const final { return os << "FormulaFunction("<<this->_formula<<")"; }
    template<class X> Void _compute(X& r, const X& x) const { r=Ariadne::evaluate(_formula,Vector<X>({x})); }
};

typedef ScalarUnivariateFormulaFunction<EffectiveNumber> EffectiveScalarUnivariateFormulaFunction;


//! A function defined by a formula
template<class Y>
struct VectorUnivariateFormulaFunction
    : VectorUnivariateFunctionMixin<VectorUnivariateFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    Vector<Formula<Y>> _formulae;

    VectorUnivariateFormulaFunction(const List< Formula<Y> >& f) : _formulae(f) { }
    VectorUnivariateFormulaFunction(const Vector< Formula<Y> >& f) : _formulae(f) { }

    ScalarUnivariateFormulaFunction<Y> operator[](SizeType i) const { return ScalarUnivariateFormulaFunction(_formulae[i]); }

    virtual SizeType result_size() const { return this->_formulae.size(); }
    virtual SizeOne argument_size() const { return SizeOne(); }
    virtual ScalarUnivariateFormulaFunction<Y>* _get(SizeType i) const { return new ScalarUnivariateFormulaFunction<Y>(this->_formulae[i]); }
    virtual VectorUnivariateFormulaFunction<Y>* _derivative(SizeOne k) const {
        return new VectorUnivariateFormulaFunction<Y>(Vector<Formula<Y>>(this->_formulae.size(),[&](SizeType i){return derivative(this->_formulae[i],k);})); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_formulae; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "VectorUnivariateFormulaFunction("<<this->result_size()<<","<<this->_formulae<<")"; }
    template<class X> Void _compute(Vector<X>& r, const X& x) const { r=Ariadne::evaluate(this->_formulae,Vector<X>({x})); }
};

typedef VectorUnivariateFormulaFunction<EffectiveNumber> EffectiveVectorUnivariateFormulaFunction;

//! A function defined by a formula
template<class Y>
struct ScalarFormulaFunction
    : ScalarMultivariateFunctionMixin<ScalarFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    SizeType _argument_size;
    Formula<Y> _formula;

    ScalarFormulaFunction(SizeType as, const Formula<Y>& f) : _argument_size(as), _formula(f) { }
    operator Formula<Y>() const { return _formula; }

    virtual SizeType argument_size() const final { return _argument_size; }
    virtual SizeOne result_size() const final { return SizeOne(); }
    virtual ScalarMultivariateFunctionInterface<P>* _derivative(SizeType j) const final { return new ScalarFormulaFunction<Y>(_argument_size,Ariadne::derivative(_formula,j)); }
    virtual OutputStream& _write(OutputStream& os) const final { return os << this->_formula; }
    virtual OutputStream& repr(OutputStream& os) const final { return os << "FormulaFunction("<<this->_argument_size<<","<<this->_formula<<")"; }
    template<class X> Void _compute(X& r, const Vector<X>& x) const { r=Ariadne::evaluate(_formula,x); }
};

typedef ScalarFormulaFunction<EffectiveNumber> EffectiveScalarFormulaFunction;

//! A vector function defined by formulae
template<class Y>
struct VectorFormulaFunction
    : VectorMultivariateFunctionMixin<VectorFormulaFunction<Y>,InformationTag<Y>>
{
    SizeType _argument_size;
    Vector< Formula<Y> > _formulae;

    VectorFormulaFunction(SizeType as, const List< Formula<Y> >& f) : _argument_size(as), _formulae(f) { }
    VectorFormulaFunction(SizeType as, const Vector< Formula<Y> >& f) : _argument_size(as), _formulae(f) { }

    ScalarFormulaFunction<Y> operator[](SizeType i) const { return ScalarFormulaFunction(_argument_size,_formulae[i]); }

    virtual SizeType result_size() const { return this->_formulae.size(); }
    virtual SizeType argument_size() const { return this->_argument_size; }
    virtual ScalarFormulaFunction<Y>* _get(SizeType i) const { return new ScalarFormulaFunction<Y>(this->_argument_size,this->_formulae[i]); }
    virtual VectorFormulaFunction<Y>* _derivative(SizeType k) const {
        return new VectorFormulaFunction<Y>(this->_argument_size, Vector<Formula<Y>>(this->_formulae.size(),[&](SizeType i){return derivative(this->_formulae[i],k);})); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_formulae; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "VectorFormulaFunction("<<this->result_size()<<","<<this->argument_size()<<","<<this->_formulae<<")"; }
    template<class X> Void _compute(Vector<X>& r, const Vector<X>& x) const { r=Ariadne::evaluate(this->_formulae,x); }
};

typedef VectorFormulaFunction<EffectiveNumber> EffectiveVectorFormulaFunction;
typedef Pair<Nat,EffectiveFormula> CoordinateFormulaPair;
typedef List<CoordinateFormulaPair> CoordinateFormulaPairs;

class NotFormulaFunctionException : public std::runtime_error {
  public:
    NotFormulaFunctionException(const String& str) : std::runtime_error(str) { }
};

//! \brief Returns \a true if the function \a f is syntactically constant in the indices \a is.
template<class Y> Bool is_constant_in(const ScalarFormulaFunction<Y>& f, const Set<Nat>& is) { return is_constant_in(f._formula,is); }
//! \brief Returns \a true if the function \a f is syntactically affine in the indices \a is.
template<class Y> Bool is_affine_in(const ScalarFormulaFunction<Y>& f, const Set<Nat>& is) { return is_affine_in(f._formula,is); }
//! \brief Returns \a true if the vector function \a f is syntactically affine in the indices \a is.
template<class Y> Bool is_affine_in(const VectorFormulaFunction<Y>& f, const Set<Nat>& is) { return is_affine_in(f._formulae,is); }
//! \brief Returns \a true if the vector function \a f is syntactically additive (possibly with multipliers) in the indices \a is.
template<class Y> Bool is_additive_in(const VectorFormulaFunction<Y>& f, const Set<Nat>& is) { return is_additive_in(f._formulae,is); }

template<class Y> Bool is_affine_in(const VectorMultivariateFunction<Y>& f, const Set<Nat>& is) {
    auto ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(f.raw_pointer());
    if (ff == nullptr) ARIADNE_THROW(NotFormulaFunctionException,"is_affine_in(f,is)","Affinity checking currently available only for formula functions.");
    return is_affine_in(ff->_formulae,is);
}
template<class Y> Bool is_additive_in(const VectorMultivariateFunction<Y>& f, const Set<Nat>& is) {
    auto ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(f.raw_pointer());
    if (ff == nullptr) ARIADNE_THROW(NotFormulaFunctionException,"is_additive_in(f,is)","Additivity checking currently available only for formula functions.");
    return is_additive_in(ff->_formulae,is);
}

inline EffectiveVectorMultivariateFunction noise_independent_component(EffectiveVectorMultivariateFunction const& function, SizeType num_inputs) {

    const EffectiveVectorFormulaFunction* ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(function.raw_pointer());
    if (ff == nullptr) ARIADNE_THROW(NotFormulaFunctionException,"noise_independent_component(f,num_inputs)","Noise independent component extraction currently available only for formula functions.");

    CoordinateFormulaPairs substitutions;
    for (auto i : range(ff->result_size(),ff->result_size()+num_inputs)) {
        substitutions.append({i,EffectiveFormula::zero()});
    }

    return EffectiveVectorFormulaFunction(function.argument_size(),simplify(substitute(ff->_formulae,substitutions)));
}

inline Vector<EffectiveVectorMultivariateFunction> input_derivatives(EffectiveVectorMultivariateFunction const& function, SizeType num_inputs) {

    Vector<EffectiveVectorMultivariateFunction> result(num_inputs);

    const EffectiveVectorFormulaFunction* ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(function.raw_pointer());
    if (ff == nullptr) ARIADNE_THROW(NotFormulaFunctionException,"input_derivatives(f,num_inputs)","Input derivatives extraction currently available only for formula functions.");

    SizeType n = function.result_size();

    for (auto j : range(num_inputs)) {
        Vector<EffectiveFormula> derivative_formulae(n);
        for (auto i : range(n)) {
            derivative_formulae[i] = simplify(derivative(ff->_formulae[i],n+j));
        }
        result[j] = EffectiveVectorFormulaFunction(function.argument_size(),derivative_formulae);
    }

    return result;
}


//------------------------ Arithmetic scalar functions  -----------------------------------//


//! A constant function f(x)=c
template<class Y, class... ARGS>
struct ConstantFunction
    : ScalarFunctionMixin<ConstantFunction<Y,ARGS...>,InformationTag<Y>,ARGS...>
{
    using P=InformationTag<Y>;
    using SIG=Real(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;

    D _domain;
    Y _value;

    //ConstantFunction(SizeType as, const Y& c) : _argument_size(as), _value(c) { }
    ConstantFunction(DomainType dom, const Y& c) : _domain(dom), _value(c) { }
    operator Y() const { return _value; }

    virtual const DomainType domain() const { return _domain; }
    virtual ArgumentSizeType argument_size() const { return _domain.dimension(); }
    virtual SizeOne result_size() const { return SizeOne(); }
    virtual ScalarFunctionInterface<P,ARGS...>* _derivative(ElementIndexType<D> j) const { return new ConstantFunction<Y,ARGS...>(_domain,Y(0)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_value; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "CF[R"<<this->argument_size()<<"]("<<_value<<")"; }
    template<class X> inline Void _compute(X& r, const ElementType<D,X>& x) const {
        r=_make_constant(_value,x); }
  private:
    template<class X> inline static X _make_constant(Y y, X const& x) { return make_constant(y,x); }
    template<class X> inline static X _make_constant(Y y, Vector<X> const& vx) { return make_constant(y,vx.zero_element()); }


};

template<class P> using ConstantUnivariateFunction = ConstantFunction<P,RealScalar>;
template<class P> using ConstantMultivariateFunction = ConstantFunction<P,RealVector>;


//! A coordinate function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=x_i\f$.
template<class P, class... ARGS>
struct CoordinateFunction
    : ScalarFunctionMixin<CoordinateFunction<P,ARGS...>,P,ARGS...>
{
    typedef Number<P> Y;
    using SIG=Real(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    DomainType _domain;
    ArgumentIndexType _index;

    //CoordinateFunction(SizeType as, SizeType i) : _argument_size(as), _index(i) { }
    CoordinateFunction(DomainType dom, ArgumentIndexType i) : _domain(dom), _index(i) { }
    ArgumentIndexType index() const { return _index; }

    virtual const DomainType domain() const { return _domain; }
    virtual ArgumentSizeType argument_size() const { return _domain.dimension(); }
    virtual SizeOne result_size() const { return SizeOne(); }
    virtual ScalarFunctionInterface<P,ARGS...>* _derivative(ArgumentIndexType j) const {
        if(j==_index) { return new ConstantFunction<Y,ARGS...>(_domain,Y(1)); }
        else { return new ConstantFunction<Y,ARGS...>(_domain,Y(0)); } }
    virtual OutputStream& _write(OutputStream& os) const { return os << "x"<<this->_index; }
    virtual OutputStream& _repr(OutputStream& os) const { return os << "IF[R"<<this->argument_size()<<"](x"<<this->_index<<")"; }
    template<class X> inline Void _compute(X& r, const Vector<X>& x) const { r=x[_index]; }
    template<class X> inline Void _compute(X& r, const Scalar<X>& x) const { r=x; }
};


//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
template<class P, class... ARGS>
struct UnaryFunction
    : ScalarFunctionMixin< UnaryFunction<P,ARGS...>, P,ARGS...>
{
    typedef Number<P> Y;
    using SIG=RealVector(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    UnaryFunction(const UnaryElementaryOperator& op, const ScalarFunction<P,ARGS...>& arg)
        : _op(op), _arg(arg) { }
    virtual UnaryFunction<P,ARGS...>* clone() const { return new UnaryFunction<P,ARGS...>(*this); }
    virtual const DomainType domain() const { return this->_arg.domain(); }
    virtual ArgumentSizeType argument_size() const { return this->_arg.argument_size(); }
    virtual SizeOne result_size() const { return SizeOne(); }

    virtual ScalarFunctionInterface<P,ARGS...>* _derivative(ArgumentIndexType j) const {
        return static_cast<const ScalarFunctionInterface<P,ARGS...>&>(this->derivative(j))._clone();
    }

    virtual ScalarFunction<P,ARGS...> derivative(ArgumentIndexType j) const {
        return _op.accept([&](auto op){
            if constexpr(Same<decltype(op),Abs>) { assert(false); return _arg.derivative(j); }
            else { return op.derivative(this->_arg,_arg.derivative(j)); } } );
    }

    virtual OutputStream& _repr(OutputStream& os) const {
        return os << "UF[R" << this->argument_size() << "](" << *this << ")"; }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << _op << '(' << _arg << ')'; }

    template<class X> inline Void _compute(X& r, const ElementType<D,X>& x) const {
        r=_op(_arg.evaluate(x)); }

    UnaryElementaryOperator _op;
    ScalarFunction<P,ARGS...> _arg;
};


template<class P, class... ARGS> ScalarFunction<P,ARGS...> sqr(const ScalarFunction<P,ARGS...>& f) {
    return ScalarFunction<P,ARGS...>(new UnaryFunction<P,ARGS...>(OperatorCode::SQR,f)); }
template<class P, class... ARGS> ScalarFunction<P,ARGS...> sqrt(const ScalarFunction<P,ARGS...>& f) {
    return ScalarFunction<P,ARGS...>(new UnaryFunction<P,ARGS...>(OperatorCode::SQRT,f)); }
template<class P, class... ARGS> ScalarFunction<P,ARGS...> sin(const ScalarFunction<P,ARGS...>& f) {
    return ScalarFunction<P,ARGS...>(new UnaryFunction<P,ARGS...>(OperatorCode::SIN,f)); }
template<class P, class... ARGS> ScalarFunction<P,ARGS...> cos(const ScalarFunction<P,ARGS...>& f) {
    return ScalarFunction<P,ARGS...>(new UnaryFunction<P,ARGS...>(OperatorCode::COS,f)); }

template<class P, class... ARGS>
struct BinaryFunction
    : ScalarFunctionMixin< BinaryFunction<P,ARGS...>, P,ARGS... >
{
    typedef Number<P> Y;
    using SIG=RealVector(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    BinaryFunction(BinaryArithmeticOperator op, const ScalarFunction<P,ARGS...>& arg1, const ScalarFunction<P,ARGS...>& arg2)
        : BinaryFunction(BinaryElementaryOperator(op.code()),arg1,arg2) { }
    BinaryFunction(BinaryElementaryOperator op, const ScalarFunction<P,ARGS...>& arg1, const ScalarFunction<P,ARGS...>& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) { ARIADNE_ASSERT_MSG(arg1.argument_size()==arg2.argument_size(),"op='"<<op<<"', arg1="<<arg1<<", arg2="<<arg2); }
    virtual BinaryFunction<P,ARGS...>* clone() const { return new BinaryFunction<P,ARGS...>(*this); }
    virtual const DomainType domain() const {
        return intersection(this->_arg1.domain(),this->_arg2.domain()); }
    virtual ArgumentSizeType argument_size() const {
        return this->_arg1.argument_size(); }
    virtual SizeOne result_size() const {
        return SizeOne(); }

    virtual ScalarFunctionInterface<P,ARGS...>* _derivative(ElementIndexType<D> j) const {
        return static_cast<const ScalarFunctionInterface<P,ARGS...>&>(this->derivative(j))._clone(); }

    virtual ScalarFunction<P,ARGS...> derivative(ElementIndexType<D> j) const {
        return _op.accept([&](auto op){
            if constexpr(Same<decltype(op),Max> || Same<decltype(op),Min>) { assert(false); return _arg1.derivative(j); }
            else { return op.derivative(_arg1,_arg1.derivative(j),_arg2,_arg2.derivative(j)); } }); }

    virtual OutputStream& repr(OutputStream& os) const {
        return os << "BF[R" << this->argument_size() << "](" << *this << ")"; }
    virtual OutputStream& _write(OutputStream& os) const {
        if(_op.code()==OperatorCode::ADD || _op.code()==OperatorCode::SUB) { return os << '(' << _arg1 << symbol(_op.code()) << _arg2 << ')'; }
        else { return os << _arg1 << symbol(_op.code()) << _arg2; } }

    template<class X> inline Void _compute(X& r, const ElementType<D,X>& x) const {
        r=_op(_arg1.evaluate(x),_arg2.evaluate(x)); }

    BinaryElementaryOperator _op;
    ScalarFunction<P,ARGS...> _arg1;
    ScalarFunction<P,ARGS...> _arg2;
};


// \brief The power function \f$(x,n)\mapsto x^n\f$.
template<class P, class... ARGS>
class GradedFunction
    : public ScalarFunctionMixin< GradedFunction<P,ARGS...>, P,ARGS... >
{
    typedef Number<P> Y;
    using SIG=RealVector(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;
  public:
    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    GradedFunction(GradedElementaryOperator op, const ScalarFunction<P,ARGS...>& arg1, const Int& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {  }
    virtual GradedFunction<P,ARGS...>* clone() const { return new GradedFunction<P,ARGS...>(*this); }
    virtual const DomainType domain() const {
        return this->_arg1.domain(); }
    virtual ArgumentSizeType argument_size() const {
        return this->_arg1.argument_size(); }
    virtual SizeOne result_size() const {
        return SizeOne(); }

    virtual ScalarFunctionInterface<P,ARGS...>* _derivative(ElementIndexType<D> j) const {
        ARIADNE_NOT_IMPLEMENTED;
    }

    virtual ScalarFunction<P,ARGS...> derivative(ElementIndexType<D> j) const {
        return _op.accept([&](auto op){ return op.derivative(_arg1, _arg1.derivative(j), _arg2);});
    }

    virtual OutputStream& repr(OutputStream& os) const {
        return os << "GF["<<this->argument_size()<<"]("<< *this <<")"; }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << _op.code() << "(" << _arg1 << "," << _arg2 << ")"; }

    template<class X> inline Void _compute(X& r, const ElementType<D,X>& x) const {
        r=_op(_arg1.evaluate(x),_arg2); }

    GradedElementaryOperator _op;
    ScalarFunction<P,ARGS...> _arg1;
    Int _arg2;
};

template<class P> using UnaryMultivariateFunction = UnaryFunction<P,RealVector>;
template<class P> using BinaryMultivariateFunction = BinaryFunction<P,RealVector>;
template<class P> using GradedMultivariateFunction = GradedFunction<P,RealVector>;

//------------------------ Vector of Scalar functions  -----------------------------------//

template<class P, class... ARGS>
class NonResizableScalarFunction : public ScalarFunction<P,ARGS...> {
  public:
    NonResizableScalarFunction<P,ARGS...>& operator=(const ScalarFunction<P,ARGS...>& f) {
        ARIADNE_ASSERT_MSG(this->domain()==f.domain(), "this->domain()="<<this->domain()<<", f.domain()="<<f.domain()<<"\n\n*this="<<*this<<"\nf="<<f<<"\n");
        this->ScalarFunction<P,ARGS...>::operator=(f);
        return *this;
    }
};

template<class P, class... ARGS>
struct VectorOfScalarFunction
    : VectorFunctionMixin<VectorOfScalarFunction<P,ARGS...>,P,ARGS...>
    , public virtual VectorOfFunctionInterface<P,ARGS...>
{
    using SIG=RealVector(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;
    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    VectorOfScalarFunction(SizeType rs, SizeType as)
        : VectorOfScalarFunction(rs, ScalarFunction<P,ARGS...>(as)) { }
    VectorOfScalarFunction(SizeType rs, DomainType dom)
        : VectorOfScalarFunction(rs, ScalarFunction<P,ARGS...>(dom)) { }
    VectorOfScalarFunction(SizeType rs, const ScalarFunction<P,ARGS...>& f)
        : _dom(f.domain()), _vec(rs,f) { }
    VectorOfScalarFunction(const Vector<ScalarFunction<P,ARGS...>>& vsf)
        : _dom(vsf.zero_element().domain()), _vec(vsf) { }

    Void set(SizeType i, ScalarFunction<P,ARGS...> f) {
        if(this->argument_size()==0u) { this->_dom=f.domain(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    ScalarFunction<P,ARGS...> get(SizeType i) const {
        return this->_vec[i]; }

    virtual SizeType result_size() const final {
        return _vec.size(); }
    virtual ArgumentSizeType argument_size() const final {
        return _dom.dimension(); }
    virtual DomainType const domain() const final {
        return _dom; }

    virtual ScalarFunctionInterface<P,ARGS...>* _get(SizeType i) const final {
        return this->_vec[i].raw_pointer()->_clone(); }
    virtual Void _set(SizeType i, const ScalarFunctionInterface<P,ARGS...>* sf) final {
        this->_vec[i]=ScalarFunction<P,ARGS...>(sf->_clone()); }
    virtual VectorFunctionInterface<P,ARGS...>* _derivative(ArgumentIndexType k) const {
        return new VectorOfScalarFunction<P,ARGS...>(Vector<ScalarFunction<P,ARGS...>>(this->_vec.size(),[&](SizeType i){return derivative(this->_vec[i],k);})); }

    const ScalarFunction<P,ARGS...> operator[](SizeType i) const {
        return this->_vec[i]; }

    NonResizableScalarFunction<P,ARGS...>& operator[](SizeType i) {
        return static_cast<NonResizableScalarFunction<P,ARGS...>&>(this->_vec[i]); }

    virtual OutputStream& _write(OutputStream& os) const {
        os << "[";
        for(SizeType i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->_write(os); }
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
    Vector<ScalarFunction<P,ARGS...>> _vec;

};

template<class P> using VectorOfScalarUnivariateFunction = VectorOfScalarFunction<P,RealScalar>;
template<class P> using VectorOfScalarMultivariateFunction = VectorOfScalarFunction<P,RealVector>;

template<class P, class... ARGS>
struct FunctionElement
    : ScalarFunctionMixin<FunctionElement<P,ARGS...>,P,ARGS...>
{
    using SIG=RealScalar(ARGS...);
    using D=typename SignatureTraits<SIG>::DomainType;

    typedef D DomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    FunctionElement(const VectorFunction<P,ARGS...>& f, SizeType i)
        : _f(f), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }

    virtual ArgumentSizeType argument_size() const { return _f.argument_size(); }
    virtual DomainType domain() const { return _f.domain(); }
    virtual OutputStream& _write(OutputStream& os) const { return os<<_f<<"["<<_i<<"]"; }
    virtual ScalarFunctionInterface<P,ARGS...>* _derivative(ArgumentIndexType j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> inline Void _compute(X& r, const ElementType<D,X>& x) const {
        r=this->_f.evaluate(x)[_i]; }

    VectorFunction<P,ARGS...> _f;
    SizeType _i;
};

//------------------------ Results of functional operations  -----------------------------------//

template<class P, class D1, class D2, class D3, class C>
struct EmbeddedFunction
    : FunctionMixin<EmbeddedFunction<P,D1,D2,D3,C>,P,ElementKind<C>(ElementKind<CartesianProductType<D1,D2,D3>>)>
{
    typedef CartesianProductType<D1,D2,D3> D;

    typedef ElementKind<D2> ARG;
    typedef ElementKind<D2> ARG2;
    typedef ElementKind<C> RES;



    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementSizeType<C> ResultSizeType;

    EmbeddedFunction(ElementSizeType<D1> as1, const Function<P,RES(ARG2)>& f2, ElementSizeType<D3> as3)
        : _dom1((as1)), _f2(f2), _dom3(as3) { ARIADNE_NOT_IMPLEMENTED; }
    EmbeddedFunction(D1 dom1, const Function<P,RES(ARG2)>& f2, D3 dom3)
        : _dom1(dom1), _f2(f2), _dom3(dom3) { }
    virtual DomainType const domain() const { return product(_dom1,_f2.domain(),_dom3); }
    virtual CodomainType const codomain() const { return _f2.codomain(); }
    virtual ArgumentSizeType argument_size() const { return _dom1.dimension()+_f2.argument_size()+_dom3.dimension(); }
    virtual ResultSizeType result_size() const { return _f2.result_size(); }
    virtual FunctionInterface<P,RES(ARG)>* _derivative(ArgumentSizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << "EmbeddedFunction( dom1="<<_dom1<<", f2="<<_f2<<", dom3="<<_dom3<<" )"; }

    template<class X> inline Void _compute(ElementType<C,X>& r, const Vector<X>& x) const {
        Vector<X> px=project(x,Range(_dom1.dimension(),_dom1.dimension()+_f2.argument_size())); r=_f2.evaluate(px); }

    D1 _dom1;
    Function<P,RES(ARG2)> _f2;
    D3 _dom3;
};


template<class P, class R, class T, class... AS>
struct ComposedFunction;

template<class P, class T, class... AS>
struct ComposedFunction<P,RealScalar,T,AS...>
    : FunctionMixin<ComposedFunction<P,RealScalar,T,AS...>,P,RealScalar(AS...)>
{
    using R=RealScalar;

    using D=DomainOfType<AS...>;
    using E=DomainOfType<T>;
    using C=DomainOfType<R>;

    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementSizeType<C> ResultSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    ComposedFunction(const Function<P,R(T)>& f, const Function<P,T(AS...)>& g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual DomainType const domain() const { return _g.domain(); }
    virtual CodomainType const codomain() const { return _f.codomain(); }
    virtual ArgumentSizeType argument_size() const { return _g.argument_size(); }
    virtual ResultSizeType result_size() const { return _f.result_size(); }
    virtual FunctionInterface<P,R(AS...)>* _derivative(ArgumentIndexType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << "ComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class X> inline Void _compute(ElementType<C,X>& r, const ElementType<D,X>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    Function<P,R(T)> _f;
    Function<P,T(AS...)> _g;
};


template<class P, class T, class... AS>
struct ComposedFunction<P,RealVector,T,AS...>
    : VectorFunctionMixin<ComposedFunction<P,RealVector,T,AS...>,P,AS...>
{
    using R=RealVector;

    using D=DomainOfType<AS...>;
    using E=DomainOfType<T>;
    using C=DomainOfType<R>;

    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementSizeType<C> ResultSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    ComposedFunction(const Function<P,R(T)>& f, const Function<P,T(AS...)>& g)
        : _f(f), _g(g) { ARIADNE_ASSERT(f.argument_size()==g.result_size()); }
    virtual DomainType const domain() const { return _g.domain(); }
    virtual CodomainType const codomain() const { return _f.codomain(); }
    virtual ArgumentSizeType argument_size() const { return _g.argument_size(); }
    virtual ResultSizeType result_size() const { return _f.result_size(); }
    virtual FunctionInterface<P,R(AS...)>* _derivative(ArgumentIndexType j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const { return os << "ComposedFunction( f="<<_f<<", g="<<_g<<" )"; }

    template<class X> inline Void _compute(ElementType<C,X>& r, const ElementType<D,X>& x) const {
        r=_f.evaluate(_g.evaluate(x)); }

    ScalarFunction<P,AS...> operator[](SizeType i) const { return compose(_f[i],_g); }

    Function<P,R(T)> _f;
    Function<P,T(AS...)> _g;
};

template<class P, class D, class C1, class C2>
struct JoinedFunction
    : VectorFunctionMixin<JoinedFunction<P,D,C1,C2>,P,ElementKind<D>>
{
    using ARG=ElementKind<D>;
    using RES1=ElementKind<C2>;
    using RES2=ElementKind<C1>;

    static_assert(Same<typename VectorFunctionMixin<JoinedFunction<P,D,C1,C2>,P,ARG>::CodomainType,CartesianProductType<C1,C2>>);

    typedef CartesianProductType<C1,C2> C;

    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementSizeType<C> ResultSizeType;
    typedef ElementIndexType<D> ArgumentIndexType;

    JoinedFunction(Function<P,RES1(ARG)> f1, Function<P,RES2(ARG)> f2)
        : _f1(f1), _f2(f2) { ARIADNE_ASSERT(f1.argument_size()==f2.argument_size()); }
    ScalarFunction<P,ARG> operator[](SizeType i) const {
        return (i<_f1.result_size()) ? _f1[i] : _f2[i-_f1.result_size()]; }

    virtual DomainType const domain() const { return intersection(_f1.domain(),_f2.domain()); }
    virtual CodomainType const codomain() const { return product(_f1.codomain(),_f2.codomain()); }
    virtual SizeType result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual ArgumentSizeType argument_size() const { return _f1.argument_size(); }
    virtual OutputStream& _write(OutputStream& os) const { return os << "JoinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }
    virtual ScalarFunctionInterface<P,ARG>* _get(SizeType i) const {
        return (i<_f1.result_size()) ? dynamic_cast<VectorOfFunctionInterface<P,ARG>const*>(_f1.raw_pointer())->_get(i)
                                     : dynamic_cast<VectorOfFunctionInterface<P,ARG>const*>(_f2.raw_pointer())->_get(i-_f1.result_size()); }
    virtual VectorFunctionInterface<P,ARG>* _derivative(ArgumentIndexType j) const { ARIADNE_NOT_IMPLEMENTED; }
    template<class X> inline Void _compute(ElementType<C,X>& r, const ElementType<D,X>& x) const {
        r=join(_f1.evaluate(x),_f2.evaluate(x)); }

    Function<P,RES1(ARG)> _f1;
    Function<P,RES2(ARG)> _f2;
};


template<class P, class D1, class D2, class C1, class C2>
class CombinedFunction
    : VectorFunctionMixin<CombinedFunction<P,D1,D2,C1,C2>,P,CartesianProductType<D1,D2>>
{
    template<class PP, class DD, class CC> using FunctionType = Function<PP,ElementKind<CC>(ElementKind<DD>)>;
    template<class PP, class DD, class CC> using FunctionInterfaceType = FunctionInterface<PP,ElementKind<CC>(ElementKind<DD>)>;

    typedef CartesianProductType<D1,D2> D;
    typedef CartesianProductType<C1,C2> C;

    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<D> ArgumentSizeType;
    typedef ElementSizeType<C> ResultSizeType;

    CombinedFunction(FunctionType<P,D1,C1> f1, FunctionType<P,D2,C2> f2)
        : _f1(f1), _f2(f2) { }
    virtual DomainType domain() const { return product(_f1.domain(),_f2.domain()); }
    virtual CodomainType codomain() const { return product(_f1.codomain(),_f2.codomain()); }
    virtual ResultSizeType result_size() const { return _f1.result_size()+_f2.result_size(); }
    virtual ArgumentSizeType argument_size() const { return _f1.argument_size()+_f2.argument_size(); }
    virtual OutputStream& _write(OutputStream& os) const { return os << "CombinedFunction( f1="<<_f1<<", f2="<<_f2<<" )"; }

    virtual FunctionInterfaceType<P,D,C>* _derivative(ElementIndexType<D> j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> inline Void _compute(ElementType<C,X>& r, const ElementType<D,X>& x) const {
        return r=join(_f1.evaluate(project(x,range(0,_f1.argument_size()))),
                      _f2.evaluate(project(x,range(_f1.argument_size(),this->argument_size())))); }

    FunctionType<P,D1,C1> _f1;
    FunctionType<P,D2,C2> _f2;
};


template<class P, class D>
class ProjectedFunction
    : VectorFunctionMixin<ProjectedFunction<P,D>,P,D>
{
    typedef ElementSizeType<D> ArgumentSizeType;

    ProjectedFunction(VectorFunction<P,D> f, Projection prj)
        : _f(f), _prj(prj) { ARIADNE_PRECONDITION(f.result_size()==prj.argument_size()); }
    virtual SizeType result_size() const { return _prj.result_size(); }
    virtual ArgumentSizeType argument_size() const { return _f.argument_size(); }
    virtual OutputStream& _write(OutputStream& os) const { return os << "ProjectedFunction( f="<<_f<<", prj="<<_prj<<" )"; }

    virtual VectorFunctionInterface<P,D>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> inline Void _compute(Vector<X>& r, const ElementType<D,X>& x) const {
        return r=_prj(_f(x)); }
    VectorFunction<P,D> _f;
    Projection _prj;
};


// A Lie deriviative \f$\nabla g\cdot f\f$.
template<class P>
struct LieDerivativeFunction
    : ScalarMultivariateFunctionMixin<LieDerivativeFunction<P>,P>
{
    //! \brief Construct the identity function in dimension \a n.
    LieDerivativeFunction(const ScalarMultivariateFunction<P>& g, const VectorMultivariateFunction<P>& f) {
        ARIADNE_ASSERT(g.argument_size()==f.argument_size());
        ARIADNE_ASSERT(f.result_size()==f.argument_size());
        _g=g; for(SizeType j=0; j!=g.argument_size(); ++j) { _dg[j]=g.derivative(j); } _f=f; }
    SizeType argument_size() const { return _g.argument_size(); }
    virtual ScalarMultivariateFunctionInterface<P>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream& os) const { return os << "LieDerivative( g="<<_g<<", f="<<_f<<" )"; }

    template<class X> inline Void _compute(X& r, const Vector<X>& x) const {
        //const Vector<R> fx=_f.evaluate(x); r=0; for(SizeType i=0; i!=_dg.size(); ++i) { r+=fx[i]+_dg[i].evaluate(x); } }
        Vector<X> fx=_f.evaluate(x);
        r=0;
        for(SizeType i=0; i!=_dg.size(); ++i) {
            r+=fx[i]+_dg[i].evaluate(x);
        }
    }

    ScalarMultivariateFunction<P> _g;
    Vector< ScalarMultivariateFunction<P> > _dg;
    VectorMultivariateFunction<P> _f;
};




} // namespace Ariadne

#endif
