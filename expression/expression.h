/***************************************************************************
 *            expression.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file expression.h
 *  \brief Internal expressions
 */

#ifndef ARIADNE_EXPRESSION_H
#define ARIADNE_EXPRESSION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/macros.h"
#include "utility/declarations.h"
#include "utility/pointer.h"
#include "utility/container.h"

#include "numeric/logical.decl.h"
#include "numeric/number.decl.h"

#include "numeric/operators.h"
#include "expression/valuation.h"
#include "expression/variables.h"
//#include "expression/operations.h"

namespace Ariadne {

class String;

template<class T> class Set;

class Identifier;

template<class T> class Variable;
template<class T> class Space;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

template<class T, class X> class Valuation;
typedef Valuation<String,String> StringValuation;
typedef Valuation<Integer,Integer> IntegerValuation;
class DiscreteValuation;
template<class X> class ContinuousValuation;

template<class X> class Vector;
template<class X> class Formula;
template<class X> class Algebra;


typedef Expression<Boolean> DiscretePredicate;
typedef Expression<Kleenean> ContinuousPredicate;
typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;

template<class T> struct DeclareExpressionOperations;


template<class X> struct ExpressionNode;

//! \ingroup ExpressionModule
//! \brief A simple expression in named variables.
//!
//! %Ariadne supports expressions of type Boolean, Kleenean, String, Integer and Real.
//! The class Real is a dummy type which can be implemented in many different ways.
//!
//! The independent variables are given string names, rather than an integer index.
//! Formulae in different variables may be combined; the variables of the resulting formula
//! are all variables occuring in all formulae.
//! Formulae may be manipulated symbolically.
//!
//! \sa Variable \sa Assignment
template<class T>
class Expression
    : public DeclareExpressionOperations<T>
{
    typedef SharedPointer<const ExpressionNode<T>> Pointer;
  public:
    typedef Real NumericType;
    typedef T ValueType;
    typedef Constant<T> ConstantType;
    typedef Variable<T> VariableType;
  public:
    explicit Expression(SharedPointer<const ExpressionNode<T>> const& eptr) : _root(eptr) { }
  public:
    //! \brief Default expression is a constant with value \c 0.
    Expression();
    //! \brief Construct an expression from a numerical value.
    Expression(const T& c);
    //! \brief Construct an expression from a named constant.
    Expression(const Constant<T>& c);
    //! \brief Construct an expression from a variable.
    Expression(const Variable<T>& v);
    //! \brief Construct a constant expression from a value.
    static Expression<T> constant(const ValueType& c);
    //! \brief Construct an expression from a name.
    static Expression<T> variable(const Identifier& c);

    //! \brief Create the zero element.
    Expression<T> create_zero() const { return Expression<T>::constant(T()); }
  public:
    const Operator& op() const;
    OperatorCode code() const;
    OperatorKind kind() const;
    const ValueType& val() const;
    const Identifier& var() const;
    const Expression<T>& arg() const;
    const Int& num() const;
    const Expression<T>& arg1() const;
    const Expression<T>& arg2() const;
    template<class A> const Expression<A>& cmp1(A* dummy=0) const;
    template<class A> const Expression<A>& cmp2(A* dummy=0) const;
    friend OutputStream& operator<<(OutputStream& os, Expression<T> const& e) { return e._write(os); }
  public:
    template<class X, EnableIf<And<IsSame<T,Real>,IsSame<X,EffectiveNumericType>>> =dummy> operator Algebra<X>() const;
  public:
    //! \brief The variables needed to compute the expression.
    Set<UntypedVariable> arguments() const;
  public:
    SharedPointer<const ExpressionNode<T>> node_ptr() const { return _root; }
    const ExpressionNode<T>* node_raw_ptr() const { return _root.operator->(); }
  private:
    OutputStream& _write(OutputStream& os) const;
  private:
    SharedPointer<const ExpressionNode<T>> _root;
};


//@{
//! \name Evaluation and related operations.
//! \related Expression

Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& q);
String evaluate(const Expression<String>& e, const StringValuation& q);
Integer evaluate(const Expression<Integer>& e, const IntegerValuation& q);
Real evaluate(const Expression<Real>& e, const ContinuousValuation<Real>& q);
Kleenean evaluate(const Expression<Kleenean>& e, const ContinuousValuation<Real>& q);

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class A> typename Logic<A>::Type evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x);

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class T> T evaluate(const Expression<T>& e, const Map<Identifier,T>& x);

//! \brief Extract the arguments of expression \a e.
template<class T> Set<Identifier> arguments(const Expression<T>& e);

//! \brief Returns \a true if the expression\a e is syntactically equal to the constant \a c.
template<class T> Bool is_constant(const Expression<T>& e, const typename Expression<T>::ValueType& c);

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable with name \a vn.
template<class T> Bool is_variable(const Expression<T>& e, const Identifier& vn);

//! \brief Returns \a true if the expression \a e is syntactically equal to the variable \a v.
template<class T> Bool is_variable(const Expression<T>& e, const Variable<T>& v);

//! \brief Simplify the expression \a e.
template<class T> Expression<T> simplify(const Expression<T>& e);

//! \brief Tests whether two expressions are identical.
template<class T> Bool identical(const Expression<T>& e1, const Expression<T>& e2);

//! \brief Returns true if the expressions are mutual negations.
//!
//! Currently can only test for pairs of the form (a1<=a2; a1>=a2),  (a1<=a2; a2<=a1)
//! or (a1>=a2; a2>=a1).
Bool opposite(Expression<Kleenean> p, Expression<Kleenean> q);

//! \brief Given \a sign when the predicate \a p is true.
Expression<Real> indicator(Expression<Kleenean> p, Sign sign=POSITIVE);

//! \brief Substitute all occurrences of variable \a v of type \c Y with constant value \a c.
template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const Variable<Y>& v, const Y& c);

//! \brief Substitute all occurrences of variable \a v of type \c Y with expression value \a se.
template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);
template<class T, class Y> Vector<Expression<T>> substitute(const Vector<Expression<T>>& e, const List< Assignment< Variable<Y>,Expression<Y> > >& a);

//! \brief Simplify the expression \a e.
Expression<Real> derivative(const Expression<Real>& e, Variable<Real> v);


//! \brief Make a formula in terms of numbered coordinates from an expression in named variables.
Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Space<Real>& spc);
Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& e, const Space<Real> spc);

//! \brief Make a function on the real line given an expression in a single argument variable.
ScalarUnivariateFunction<EffectiveTag> make_function(const Variable<Real>& v, const Expression<Real>& e);
//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables.
ScalarFunction<EffectiveTag> make_function(const Space<Real>& s, const Expression<Real>& e);
//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables.
VectorFunction<EffectiveTag> make_function(const Space<Real>& s, const Vector<Expression<Real>>& e);

//! \brief Make a function on a Euclidean domain given an ordered list including all argument variables. // DEPRECATED
ScalarFunction<EffectiveTag> make_function(const Expression<Real>& e, const Space<Real>& s);

//@}




} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_H */
