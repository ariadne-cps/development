/***************************************************************************
 *            expression.tpl.hpp
 *
 *  Copyright  2009-20  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOperatorCode::OR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "utility/standard.hpp"

#include "numeric/operators.tpl.hpp"

#include "algebra/algebra.hpp"
#include "algebra/algebra_wrapper.hpp"

#include "symbolic/constant.hpp"
#include "symbolic/variable.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/space.hpp"
#include "symbolic/valuation.hpp"

#include "function/formula.hpp"

namespace Ariadne {

using Eq = Equal;
using Neq = Unequal;

template<class SIG> struct OperatorTypedef;
template<class SIG> using OperatorType = typename OperatorTypedef<SIG>::Type;

struct CatOp {
    static constexpr OperatorCode code() { return OperatorCode::ADD; }
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    String operator()(String const& s1, String const& s2) const { return s1+s2; }
    friend OutputStream& operator<<(OutputStream& os, CatOp op) { return os << "cat"; }
};
inline String pos(String s) { return s; }


template<> struct OperatorTypedef<String(String)> { typedef OperatorVariant<> Type; };

template<> struct OperatorTypedef<Integer(Integer)> { typedef OperatorVariant<Nul,Pos,Neg,Sqr,Abs> Type; };
template<> struct OperatorTypedef<Integer(Integer,Integer)> { typedef OperatorVariant<Add,Sub,Mul,Max,Min> Type; };

template<> struct OperatorTypedef<Real(Real)> { typedef UnaryElementaryOperator Type; };
template<> struct OperatorTypedef<Real(Real,Real)> { typedef BinaryElementaryOperator Type; };
template<> struct OperatorTypedef<Real(Real,Integer)> { typedef GradedElementaryOperator Type; };
template<> struct OperatorTypedef<Real(Real,Int)> { typedef GradedElementaryOperator Type; };

template<> struct OperatorTypedef<Boolean(String,String)> { typedef OperatorVariant<Eq,Neq> Type; };
template<> struct OperatorTypedef<Boolean(Integer,Integer)> { typedef BinaryComparisonOperator Type; };
template<> struct OperatorTypedef<Kleenean(Real)> { typedef OperatorVariant<Sgn> Type; };
//template<> struct OperatorTypedef<Kleenean(Real,Real)> { typedef OperatorVariant<Less,Gtr,Leq,Geq> Type; };
template<> struct OperatorTypedef<Kleenean(Real,Real)> { typedef BinaryComparisonOperator Type; };

template<> struct OperatorTypedef<Boolean(Boolean)> { typedef UnaryLogicalOperator Type; };
template<> struct OperatorTypedef<Boolean(Boolean,Boolean)> { typedef BinaryLogicalOperator Type; };
template<> struct OperatorTypedef<Kleenean(Kleenean)> { typedef UnaryLogicalOperator Type; };
template<> struct OperatorTypedef<Kleenean(Kleenean,Kleenean)> { typedef BinaryLogicalOperator Type; };

template<class SIG> struct SymbolicTypedef;
template<class SIG> using SymbolicType = typename SymbolicTypedef<SIG>::Type;

template<class R, class A> struct SymbolicTypedef<R(A)> {
    typedef Symbolic<OperatorType<R(A)>,Expression<A>> Type; };
template<class R, class A1, class A2> struct SymbolicTypedef<R(A1,A2)> {
    typedef Symbolic<OperatorType<R(A1,A2)>,Expression<A1>,Expression<A2>> Type; };
template<class R, class... AS> struct SymbolicTypedef<R(AS...)> {
    typedef Symbolic<OperatorType<R(AS...)>,Expression<AS>...> Type; };


template<class T> struct ExpressionVariantTypedef;
template<class T> using ExpressionVariantType = typename ExpressionVariantTypedef<T>::Type;

template<> struct ExpressionVariantTypedef<String> {
    typedef Variant< Constant<String>, Variable<String> > Type;
};

template<> struct ExpressionVariantTypedef<Integer> {
    typedef Variant< Constant<Integer>, Variable<Integer>,
                     SymbolicType<Integer(Integer)>, SymbolicType<Integer(Integer,Integer)> > Type;
};

template<> struct ExpressionVariantTypedef<Real> {
    typedef Variant< Constant<Real>, Variable<Real>,
                     SymbolicType<Real(Real)>, SymbolicType<Real(Real,Real)>,
                     Symbolic<OperatorType<Real(Real,Integer)>,Expression<Real>,Int> > Type;
};

template<> struct ExpressionVariantTypedef<Boolean> {
    typedef Variant< Constant<Boolean>, Variable<Boolean>,
                     SymbolicType<Boolean(Boolean)>, SymbolicType<Boolean(Boolean,Boolean)>,
                     SymbolicType<Boolean(String,String)>, SymbolicType<Boolean(Integer,Integer)> > Type;
};

template<> struct ExpressionVariantTypedef<Kleenean> {
    typedef Variant< Constant<Kleenean>, Variable<Kleenean>,
                     SymbolicType<Kleenean(Kleenean)>, SymbolicType<Kleenean(Kleenean,Kleenean)>,
                     SymbolicType<Kleenean(Real)>, SymbolicType<Kleenean(Real,Real)> > Type;
};

namespace {
template<class T> inline Cnst _op_impl(Constant<T> const&) { Cnst op; return op; }
template<class T> inline Var _op_impl(Variable<T> const&) { return Var(); }
template<class OP, class... AS> inline OP _op_impl(Symbolic<OP,AS...> const& s) { return s._op; }
}

template<class T, class... TS> struct IndexOf<T,Variant<TS...>> { static const SizeType N=IndexOf<T,TS...>::N; };

template<class T, class VAR> constexpr decltype(auto) variant_index_of() { return IntegralConstant<SizeType,IndexOf<T,VAR>::N>(); }

template<class T> struct ExpressionNode : public ExpressionVariantType<T> {
  public:
    template<class... AS> requires Constructible<ExpressionVariantType<T>,AS...>
        ExpressionNode(AS... as) : ExpressionVariantType<T>(as...) { }

    ExpressionVariantType<T> const& base() const { return *this; }

    template<class VIS> decltype(auto) accept(VIS&& vis) const {
        return std::visit(std::forward<VIS>(vis),static_cast<ExpressionVariantType<T>const&>(*this)); }

    template<class EN> static constexpr decltype(auto) index_of() {
        return variant_index_of<EN,ExpressionVariantType<T>>(); }
    Operator op() const { return this->accept([](auto s){return Operator(_op_impl(s));}); }
};

template<class T> inline OutputStream& operator<<(OutputStream& os, const ExpressionNode<T>* e) {
    return os << static_cast<Void const*>(e);
}

template<class T> using ConstantExpressionNode = Constant<T>;
template<class T> using NamedConstantExpressionNode = ConstantExpressionNode<T>;
template<class T> using VariableExpressionNode = Variable<T>;

template<class R, class A=R> using UnaryExpressionNode = Symbolic<OperatorType<R(A)>,Expression<A>>;
template<class R, class A1=R, class A2=A1> using BinaryExpressionNode = Symbolic<OperatorType<R(A1,A2)>,Expression<A1>,Expression<A2>>;
template<class R, class A=R, class N=Int> using GradedExpressionNode = Symbolic<OperatorType<R(A,N)>,Expression<A>,N>;


template<class T> Expression<T>::Expression() : Expression(T()) { }
template<class T> Expression<T>::Expression(const T& c) : Expression(Constant<T>(c)) { }
template<class T> Expression<T>::Expression(const Constant<T>& c): _root(new ExpressionNode<T>(c)) { }
template<class T> Expression<T>::Expression(const Variable<T>& v) : _root(new ExpressionNode<T>(v)) { }

template<> Expression<String>::Expression(const String& c): _root(new ExpressionNode<String>(Constant<String>(c))) { }

template<class T> Expression<T> Expression<T>::constant(const T& c) {
    return Expression<T>(c); }
template<class T> Expression<T> Expression<T>::constant(const Constant<T>& c) {
    return Expression<T>(c); }
template<class T> Expression<T> Expression<T>::variable(const Identifier& v) {
    return Expression<T>(Variable<T>(v)); }

template<class T> Operator Expression<T>::op() const {
    return this->node_ptr()->op(); }
template<class T> OperatorCode Expression<T>::code() const {
    return node_ptr()->op().code(); }
template<class T> OperatorKind Expression<T>::kind() const {
    return node_ptr()->op().kind(); }
template<class T> const T& Expression<T>::val() const {
    return std::get<Constant<T>>(node_ref()); }
template<class T> const Identifier& Expression<T>::var() const {
    return std::get<VariableExpressionNode<T>>(node_ref()).name(); }
template<class T> const Expression<T>& Expression<T>::arg() const {
    if constexpr (IsSame<T,Real>::value) { if (auto* gn = std::get_if<GradedExpressionNode<T>>(&node_ref())) { return gn->_arg; } }
    if constexpr (not IsSame<T,String>::value) { return std::get<UnaryExpressionNode<T>>(node_ref())._arg; } else { abort(); } }
template<class T> const Int& Expression<T>::num() const {
    if constexpr (IsSame<T,Real>::value) { return std::get<GradedExpressionNode<T>>(node_ref())._num; } else { abort(); } }
template<class T> const Expression<T>& Expression<T>::arg1() const {
    if constexpr (not IsSame<T,String>::value) { return std::get<BinaryExpressionNode<T>>(node_ref())._arg1; } else { abort(); } }
template<class T> const Expression<T>& Expression<T>::arg2() const {
    if constexpr (not IsSame<T,String>::value) { return std::get<BinaryExpressionNode<T>>(node_ref())._arg2; } else { abort(); } }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp1(A*) const {
    return std::get<Symbolic<OperatorType<R(A,A)>,Expression<A>,Expression<A>>>(node_ref())._arg1; }
template<class R> template<class A> const Expression<A>& Expression<R>::cmp2(A*) const {
    return std::get<Symbolic<OperatorType<R(A,A)>,Expression<A>,Expression<A>>>(node_ref())._arg2; }


template<class T> Set<UntypedVariable> Expression<T>::arguments() const {
    return this->node_ref().accept([](auto s){return _arguments(s);});
}

template<class T> OutputStream& Expression<T>::_write(OutputStream& os) const {
    return os << _default_writer(*this);
}

template<class T> Writer<Expression<T>> Expression<T>::_default_writer(new PrefixExpressionWriter<T>());

template<class T> OutputStream& PrefixExpressionWriter<T>::_write(OutputStream& os, Expression<T> const& e) const {
    e.node_ref().accept([&os](auto expr){_write_impl(os,expr);}); return os;
}

template<class T> OutputStream& InfixExpressionWriter<T>::_write(OutputStream& os, Expression<T> const& e) const {
    e.node_ref().accept([&os](auto expr){_write_impl(os,expr);}); return os;
}

template<class R> inline
Expression<R> make_expression(const Constant<R>& c) {
    return Expression<R>(std::make_shared<NamedConstantExpressionNode<R>>(c)); }
template<class R> inline
Expression<R> make_expression(const R& c) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(ConstantExpressionNode<R>(Constant<R>(c)))); }
template<class R, class A> inline
Expression<R> make_expression(OperatorType<R(A)> op, const Expression<A>& e) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(UnaryExpressionNode<R,A>(op,e))); }
template<class R, class A, class N> inline
Expression<R> make_expression(OperatorType<R(A,N)> op, const Expression<A>& e, N n) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(GradedExpressionNode<R,A,N>(op,e,n))); }
template<class R, class A1, class A2> inline
Expression<R> make_expression(OperatorType<R(A1,A2)> op, const Expression<A1>& e1, Expression<A2> e2) {
    return Expression<R>(std::make_shared<ExpressionNode<R>>(BinaryExpressionNode<R,A1,A2>(op,e1,e2))); }

template<class R, class OP, class A> inline
Expression<R> make_expression(OP op, const Expression<A> e) {
    return make_expression<R>(OperatorType<R(A)>(op),e); }
template<class R, class A, class N> inline
Expression<R> make_expression(Pow op, const Expression<A>& e, N n) {
    return make_expression<R>(OperatorType<R(A,Int)>(op),e,n); }
template<class R, class OP, class A1, class A2> inline
Expression<R> make_expression(OP op, const Expression<A1>& e1, Expression<A2> e2) {
    return make_expression<R>(OperatorType<R(A1,A2)>(op),e1,e2); }



template<class A> A evaluate(const Expression<A>& e, const Map<Identifier,A>& x);

inline Integer evaluate(const Expression<Integer>& e, const Map<Identifier,String>& x) {
    return e.node_ref().accept([&x](auto en){return evaluate_as<Integer>(en,x);});
}
inline String evaluate(const Expression<String>& e, const Map<Identifier,Integer>& x) {
    return e.node_ref().accept([&x](auto en){return evaluate_as<String>(en,x);});
}


template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x) {
    typedef typename Logic<A>::Type R;
    return e.node_ref().accept([&x](auto en){return evaluate_as<R>(en,x);});
}

template<class T> T evaluate(const Expression<T>& e, const Map<Identifier,T>& x) {
    return e.node_ref().accept([&x](auto en){return evaluate_as<T>(en,x);});
}

template<class T> T evaluate(const Expression<T>& e, const Valuation<T>& x) {
    return evaluate(e,x.values());
}

template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Valuation<A>& x) {
    return evaluate(e,x.values());
}


template<class T> Set<Identifier> arguments(const Expression<T>& e) {
    Set<UntypedVariable> arg_vars = e.arguments();
    Set<Identifier> arg_names;
    for (auto var : arg_vars) { arg_names.insert(var.name()); }
    return arg_names;
}


template<class T> Bool is_constant_in(const Expression<T>& e, const Set<Variable<T>>& spc) {
    return e.node_ref().accept([&spc](auto en){return is_constant_in(en,spc);});
}


namespace {
template<class X, class Y> Expression<X> _substitute(const Constant<X>& c, const Variable<Y>& v, const Expression<Y>& s) { return c; }
template<class X, class Y> Expression<X> _substitute(const Variable<X>& var, const Variable<Y>& v, const Expression<Y>& s) { return var; }
template<class X> Expression<X> _substitute(const Variable<X>& var, const Variable<X>& v, const Expression<X>& s) {
    if (var==v) { return s; } else { return var; } }
template<class X, class Y, class OP, class E> Expression<X> _substitute(const Symbolic<OP,E>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e._op(substitute(e._arg,v,s)); }
template<class X, class Y, class OP, class E1, class E2> Expression<X> _substitute(const Symbolic<OP,E1,E2>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e._op(substitute(e._arg1,v,s),substitute(e._arg2,v,s)); }
template<class X, class Y, class OP, class E> Expression<X> _substitute(const Symbolic<OP,E,Int>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e._op(substitute(e._arg,v,s),e._num); }
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return e.node_ref().accept([&v,&s](auto en){return _substitute<X>(en,v,s);});
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c) {
    return substitute(e,v,Expression<Y>::constant(c));
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Assignment<Variable<Y>,Expression<Y>>& a) {
    return substitute(e,a.lhs,a.rhs);
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const List<Assignment<Variable<Y>,Expression<Y>>>& a) {
    Expression<X> r=e;
    for(SizeType i=0; i!=a.size(); ++i) {
        r=substitute(r,a[i]);
    }
    return r;
}

template<class X, class Y> Vector<Expression<X>> substitute(const Vector<Expression<X>>& e, const List< Assignment< Variable<Y>, Expression<Y> > >& a) {
    Vector<Expression<X>> r(e.size());
    for(SizeType i=0; i!=e.size(); ++i) {
        r[i]=substitute(e[i],a);
    }
    return r;
}



inline Bool same(Kleenean k1, Kleenean k2) { return definitely(k1==k2); }



template<class T> Bool identical(const Expression<T>& e1, const Expression<T>& e2)
{
    if(e1.node_raw_ptr()==e2.node_raw_ptr()) { return true; }
    if(e1.op()!=e2.op()) { return false; }
    return e1.node_ref().accept([&e2](auto e1n){return e2.node_ref().accept([&e1n](auto e2n){return _identical_dispatch(e1n,e2n);});});
}


namespace {

template<class X> inline Expression<X> _simplify(const Expression<X>& e) {
    return e;
}

template<class OP, class A1, class A2> decltype(auto) _simpl(OP op, Constant<A1> a1, Constant<A2> a2) {
    auto r=op(a1.value(),a2.value()); return Constant<decltype(r)>(r); }
template<class OP, class A> decltype(auto) _simpl(OP op, Constant<A> a) {
    auto r=op(a.value()); return Constant<decltype(r)>(r); }

template<class E, class T> inline E simplify_variant(const ExpressionVariantType<T>& e) {
    return e.accept([](auto en){return _simplify_node<E>(en);});
}

template<class E> inline E _simplify_node(const UnaryExpressionNode<Kleenean,Real>& e) {
    return make_expression<Kleenean>(e.op(),simplify(e.arg()));
}
template<class E> inline E _simplify_node(const BinaryExpressionNode<Kleenean,Real,Real>& e) {
    return make_expression<Kleenean>(e.op(),simplify(e.arg1()),simplify(e.arg2()));
}

Expression<Real> _simplify(const Expression<Real>& e) {
    // Need to dispatch here to allow use of specialisation for UnaryExpressionNode<T>, which correctly handles constants
    return e.node_ref().accept([](auto en){return _simplify_node<Expression<Real>>(en);});
}

Expression<Kleenean> _simplify(const Expression<Kleenean>& e) {
    return e.node_ref().accept([](auto en){return _simplify_node<Expression<Kleenean>>(en);});
}

} // namespace

template<class X> Expression<X> simplify(const Expression<X>& e) {
    return _simplify(e);
}



template<class T> Bool _before(Constant<T> const& c1, Constant<T> const& c2) {
    return decide(c1<c2); }
template<class T> Bool _before(Variable<T> const& v1, Variable<T> const& v2) {
    return v1.name() < v2.name(); }
template<class OP, class A> Bool _before(Symbolic<OP,A> const& s1, Symbolic<OP,A> const& s2) {
    if (s1.op().code() == s2.op().code()) { return before(s1.arg(),s2.arg()); }
    else { return s1.op().code() < s2.op().code(); } }
template<class OP, class A1, class A2> Bool _before(Symbolic<OP,A1,A2> const& s1, Symbolic<OP,A1,A2> const& s2) {
    if (s1.op().code() == s2.op().code()) {
        if (identical(s1.arg1(),s2.arg1())) { return before(s1.arg2(),s2.arg2()); }
        else { return before(s1.arg1(),s2.arg1()); }
    } else { return s1.op().code() < s2.op().code(); } }
template<class OP, class A> Bool _before(Symbolic<OP,A,Int> const& s1, Symbolic<OP,A,Int> const& s2) {
    if (s1.op().code() == s2.op().code()) {
        if (s1.num() == s2.num()) { return before(s1.arg(),s2.arg()); }
        else { return s1.num() < s2.num(); }
    } else { return s1.op().code() < s2.op().code(); } }

template<class T> Bool before(Expression<T> const& e1, Expression<T> const& e2) {
    if (e1.node_ref().index() == e2.node_ref().index()) {
        return e1.node_ref().accept([&e2](auto e1n){return _before(e1n,std::get<decltype(e1n)>(e2.node_ref()));});
    } else {
        return e1.node_ref().index() < e2.node_ref().index();
    }
}

struct ExpressionComparator {
    template<class A1, class A2> decltype(auto) operator() (A1&& a1, A2&& a2) const { return before(std::forward<A1>(a1),std::forward<A2>(a2)); }
};

template<class T, class CMP=ExpressionComparator> using ExpressionSet = std::set<Expression<T>,CMP>;

template<class T> class CommonSubroutineEliminator {
    ExpressionSet<T> _cache;
  public:
    CommonSubroutineEliminator() : _cache() { }
    Expression<T> eliminate_common_subexpressions(const Expression<T>& e);
    Expression<T> operator()(const ConstantExpressionNode<T>&, const Expression<T>& e) { _cache.insert(e); return e; }
    Expression<T> operator()(const VariableExpressionNode<T>&, const Expression<T>& e) { _cache.insert(e); return e; }
    Expression<T> operator()(const BinaryExpressionNode<T>& s, const Expression<T>& e) {
        auto new_arg1 = eliminate_common_subexpressions(s.arg1());
        auto new_arg2 = eliminate_common_subexpressions(s.arg2());
        if(new_arg1.node_raw_ptr() == s.arg1().node_raw_ptr() && new_arg2.node_raw_ptr() == s.arg2().node_raw_ptr()) {
            _cache.insert(e); return e;
        } else {
            Expression<T> new_e=make_expression<T>(s.op(),new_arg1,new_arg2);
            _cache.insert(new_e); return new_e;
        }
    }
    Expression<T> operator()(const UnaryExpressionNode<T>& s, const Expression<T>& e) {
        auto new_arg = eliminate_common_subexpressions(s.arg());
        if(new_arg.node_raw_ptr() == s.arg().node_raw_ptr()) {
            _cache.insert(e); return e;
        } else {
            Expression<T> new_e=make_expression<T>(s.op(),new_arg);
            _cache.insert(new_e); return new_e;
        }
    }
    Expression<T> operator()(const GradedExpressionNode<T>& s, const Expression<T>& e) {
        auto new_arg = eliminate_common_subexpressions(s.arg());
        if(new_arg.node_raw_ptr() == s.arg().node_raw_ptr()) {
            _cache.insert(e); return e;
        } else {
            Expression<T> new_e=make_expression<T>(s.op(),new_arg,s.num());
            _cache.insert(new_e); return new_e;
        }
    }
};
template<class T> Expression<T> CommonSubroutineEliminator<T>::eliminate_common_subexpressions(const Expression<T>& e) {
    auto iter=this->_cache.find(e);
    if (iter!=_cache.end()) {
        return *iter; }
    else {
        return e.node_ref().accept([&](auto en){return this->operator()(en,e);});
    }
}

template<class T> Void eliminate_common_subexpressions(Expression<T>& e)
{
    e=CommonSubroutineEliminator<T>().eliminate_common_subexpressions(e);
}

template<class T> Void eliminate_common_subexpressions(Vector<Expression<T>>& es)
{
    CommonSubroutineEliminator<T> simplifier;
    for(SizeType i=0; i!=es.size(); ++i) { es[i]=simplifier.eliminate_common_subexpressions(es[i]); }
}

struct ExpressionPtrComparator {
    template<class A1, class A2> decltype(auto) operator() (A1&& a1, A2&& a2) const { return a1.node_raw_ptr() < a2.node_raw_ptr(); }
};


template<class T, class CMP, Bool distinct> class NodeCounter {
    ExpressionSet<T,CMP> cache; SizeType count;
  public:
    NodeCounter() : cache(), count(0u) { }
    SizeType count_nodes (Expression<T> const& e);
  private: public:
    Void operator() (ConstantExpressionNode<T> const& e) { }
    Void operator() (VariableExpressionNode<T> const& e) { }
    Void operator() (BinaryExpressionNode<T> const& e) { this->count_nodes(e.arg1()); this->count_nodes(e.arg2()); }
    Void operator() (UnaryExpressionNode<T> const& e) { this->count_nodes(e.arg()); }
    Void operator() (GradedExpressionNode<T> const& e) { this->count_nodes(e.arg()); }
};

template<class T, class CMP, Bool distinct> auto NodeCounter<T,CMP,distinct>::count_nodes(Expression<T> const& e) -> SizeType {
    auto iter=cache.find(e);
    if (!distinct || iter==cache.end()) {
        cache.insert(e);
        e.node_ref().accept(*this);
        count++;
    }
    return count;
}


template<class T> SizeType count_nodes(Expression<T> const& e) {
    return NodeCounter<T,ExpressionComparator,false>().count_nodes(e);
}

template<class T> SizeType count_distinct_nodes(Expression<T> const& e) {
    return NodeCounter<T,ExpressionComparator,true>().count_nodes(e);
}

template<class T> SizeType count_distinct_node_pointers(Expression<T> const& e) {
    return NodeCounter<T,ExpressionPtrComparator,true>().count_nodes(e);
}


template<class T> Bool is_constant(const Expression<T>& e, const SelfType<T>& c) {
    auto* ce = std::get_if<Constant<T>>(&e.node_ref()); return ce && decide(ce->value()==c);
}

template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v) {
    auto* ve = std::get_if<Variable<T>>(&e.node_ref()); return ve && *ve==v;
}


} // namespace Ariadne
