/***************************************************************************
 *            function_interface.h
 *
 *  Copyright 2008  Pieter Collins
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

/*! \file function_interface.h
 *  \brief Interface for functions for which derivatives can be computed.
 */
#ifndef ARIADNE_FUNCTION_INTERFACE_H
#define ARIADNE_FUNCTION_INTERFACE_H

#include <iosfwd>

namespace Ariadne {

typedef std::ostream OutputStream;

static const int SMOOTH=255;

typedef void Void;
typedef unsigned int Nat;
typedef int Int;

class Float;
class Interval;
class Real;

class ApproximateFloat;
class ValidatedFloat;
class ExactFloat;

typedef Float ApproximateNumberType;
typedef Interval ValidatedNumberType;
typedef Real EffectiveNumberType;
typedef Float RawNumberType;

typedef ApproximateNumberType ApproximateNumber;
typedef ValidatedNumberType ValidatedNumber;
typedef EffectiveNumberType EffectiveNumber;

//typedef ApproximateFloat ApproximateNumberType;
//typedef ValidatedFloat ValidatedNumberType;
//typedef Real EffectiveNumberType;

struct ExactTag { };
typedef EffectiveNumberType EffectiveTag;
typedef ValidatedNumberType ValidatedTag;
typedef ApproximateNumberType ApproximateTag;

template<class I> struct CanonicalNumberTypedef;
template<> struct CanonicalNumberTypedef<ExactTag> { typedef EffectiveNumberType Type; };
template<> struct CanonicalNumberTypedef<EffectiveTag> { typedef EffectiveNumberType Type; };
template<> struct CanonicalNumberTypedef<ValidatedTag> { typedef ValidatedNumberType Type; };
template<> struct CanonicalNumberTypedef<ApproximateTag> { typedef ApproximateNumberType Type; };
template<class I> using CanonicalNumberType = typename CanonicalNumberTypedef<I>::Type;

class Interval;
class Box;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class TaylorModel;
template<class X> class Formula;
template<class X> class Algebra;

template<class X> class ScalarFunctionInterface;
template<class X> class VectorFunctionInterface;

typedef ScalarFunctionInterface<ApproximateTag> ApproximateScalarFunctionInterface;
typedef ScalarFunctionInterface<ValidatedTag> ValidatedScalarFunctionInterface;
typedef ScalarFunctionInterface<EffectiveTag> EffectiveScalarFunctionInterface;

typedef VectorFunctionInterface<ApproximateTag> ApproximateVectorFunctionInterface;
typedef VectorFunctionInterface<ValidatedTag> ValidatedVectorFunctionInterface;
typedef VectorFunctionInterface<EffectiveTag> EffectiveVectorFunctionInterface;

template<class X> class ScalarFunction;
typedef ScalarFunction<ApproximateTag> ApproximateScalarFunction;
typedef ScalarFunction<ValidatedTag> ValidatedScalarFunction;
typedef ScalarFunction<EffectiveTag> EffectiveScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;

template<class X> class VectorFunction;
typedef VectorFunction<ApproximateTag> ApproximateVectorFunction;
typedef VectorFunction<ValidatedTag> ValidatedVectorFunction;
typedef VectorFunction<EffectiveTag> EffectiveVectorFunction;
typedef EffectiveVectorFunction RealVectorFunction;

template<>
class ScalarFunctionInterface<Void>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef Nat SizeType;

    //! \brief Virtual destructor.
    virtual ~ScalarFunctionInterface() { };
    //! \brief The number of arguments to the expression.
    virtual SizeType argument_size() const = 0;

    //! \brief Write a full version to an output stream.
    virtual std::ostream& repr(std::ostream& os) const = 0;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
  public:
    virtual ScalarFunctionInterface<Void>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{F}^n\rightarrow\mathbb{F}\f$ which can only be evaluated approximately.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<ApproximateTag>
    : public ScalarFunctionInterface<Void>
{
  public:
    //! \brief Return a copy of the function.
    inline ScalarFunction<ApproximateTag> clone() const;
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual ApproximateNumberType evaluate(const Vector<ApproximateNumberType>& x) const = 0;
    inline ApproximateNumberType operator() (const Vector<ApproximateNumberType>& x) const;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Differential<ApproximateNumberType> evaluate(const Vector< Differential<ApproximateNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual TaylorModel<ApproximateNumberType> evaluate(const Vector< TaylorModel<ApproximateNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae to create the composed function.
    virtual Formula<ApproximateNumberType> evaluate(const Vector< Formula<ApproximateNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of elements of an algebra.
    virtual Algebra<ApproximateNumberType> evaluate(const Vector< Algebra<ApproximateNumberType> >& x) const = 0;

    Vector<ApproximateNumberType> gradient(const Vector<ApproximateNumberType>& x) const;
    Differential<ApproximateNumberType> differential(const Vector<ApproximateNumberType>& x, Nat d) const;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<ApproximateNumberType> derivative(Nat i) const;
  private:
    virtual ScalarFunctionInterface<ApproximateNumberType>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<ApproximateNumberType>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{I}^n\rightarrow\mathbb{I}\f$ which can be evaluated over intervals.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<ValidatedTag>
    : public ScalarFunctionInterface<ApproximateTag>
{
  public:

    inline ScalarFunction<ValidatedTag> clone() const;

    using ScalarFunctionInterface<ApproximateTag>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual ValidatedNumberType evaluate(const Vector<ValidatedNumberType>& x) const = 0;
    inline ValidatedNumberType operator() (const Vector<ValidatedNumberType>& x) const;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Differential<ValidatedNumberType> evaluate(const Vector< Differential<ValidatedNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual TaylorModel<ValidatedNumberType> evaluate(const Vector< TaylorModel<ValidatedNumberType> >& x) const = 0;

    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<ValidatedNumberType> evaluate(const Vector< Formula<ValidatedNumberType> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Algebra<ValidatedNumberType> evaluate(const Vector< Algebra<ValidatedNumberType> >& x) const = 0;

    using ScalarFunctionInterface<ApproximateTag>::gradient;
    Vector<ValidatedNumberType> gradient(const Vector<ValidatedNumberType>& x) const;
    using ScalarFunctionInterface<ApproximateTag>::differential;
    Differential<ValidatedNumberType> differential(const Vector<ValidatedNumberType>& x, Nat d) const;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<ValidatedTag> derivative(Nat i) const;
  private:
    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunctionInterface<ValidatedTag>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<ValidatedTag>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<EffectiveTag>
    : public ScalarFunctionInterface<ValidatedTag>
{
  public:
    inline ScalarFunction<EffectiveTag> clone() const;

    using ScalarFunctionInterface<ValidatedTag>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual EffectiveNumberType evaluate(const Vector<EffectiveNumberType>& x) const = 0;
    inline EffectiveNumberType operator() (const Vector<EffectiveNumberType>& x) const;
    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<EffectiveNumberType> evaluate(const Vector< Formula<EffectiveNumberType> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Algebra<EffectiveNumberType> evaluate(const Vector< Algebra<EffectiveNumberType> >& x) const = 0;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<EffectiveTag> derivative(Nat i) const;
  private:
    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunctionInterface<EffectiveTag>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<EffectiveTag>* _clone() const = 0;
};

//! \relates ScalarFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const ScalarFunctionInterface<Void>& f) {
    return f.write(os);
}


//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Void>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef Nat SizeType;

    //! \brief Virtual destructor.
    virtual ~VectorFunctionInterface() { };

    //! \brief The number of arguments to the function.
    virtual SizeType argument_size() const = 0;
    //! \brief The number of result variables of the function.
    virtual SizeType result_size() const = 0;

    //! \brief Write a full version to an output stream.
    virtual std::ostream& repr(std::ostream& os) const = 0;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
  public:
    virtual VectorFunctionInterface<Void>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<ApproximateTag>
    : public VectorFunctionInterface<Void>
{
  public:
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Vector<ApproximateNumberType> evaluate(const Vector<ApproximateNumberType>& x) const = 0;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Vector< Differential<ApproximateNumberType> > evaluate(const Vector< Differential<ApproximateNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual Vector< TaylorModel<ApproximateNumberType> > evaluate(const Vector< TaylorModel<ApproximateNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae to create the composed function.
    virtual Vector< Formula<ApproximateNumberType> > evaluate(const Vector< Formula<ApproximateNumberType> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<ApproximateNumberType> > evaluate(const Vector< Algebra<ApproximateNumberType> >& x) const = 0;

    Matrix<ApproximateNumberType> jacobian(const Vector<ApproximateNumberType>& x) const;
    Vector< Differential<ApproximateNumberType> > differentials(const Vector<ApproximateNumberType>& x, Nat d) const;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<ApproximateTag> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<ApproximateTag>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<ApproximateTag>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\I^n\rightarrow\I^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<ValidatedTag>
    : public VectorFunctionInterface<ApproximateTag>
{
  public:
    using VectorFunctionInterface<ApproximateTag>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Vector<ValidatedNumberType> evaluate(const Vector<ValidatedNumberType>& x) const = 0;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Vector< Differential<ValidatedNumberType> > evaluate(const Vector< Differential<ValidatedNumberType> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual Vector< TaylorModel<ValidatedNumberType> > evaluate(const Vector< TaylorModel<ValidatedNumberType> >& x) const = 0;

    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<ValidatedNumberType> > evaluate(const Vector< Formula<ValidatedNumberType> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<ValidatedNumberType> > evaluate(const Vector< Algebra<ValidatedNumberType> >& x) const = 0;

    using VectorFunctionInterface<ApproximateTag>::jacobian;
    Matrix<ValidatedNumberType> jacobian(const Vector<ValidatedNumberType>& x) const;
    using VectorFunctionInterface<ApproximateTag>::differentials;
    Vector< Differential<ValidatedNumberType> > differentials(const Vector<ValidatedNumberType>& x, Nat d) const;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<ValidatedTag> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<ValidatedTag>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<ValidatedTag>* _clone() const = 0;

};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\R^n\rightarrow\R^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<EffectiveTag>
    : public VectorFunctionInterface<ValidatedTag>
{
  public:
    using VectorFunctionInterface<ValidatedTag>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual Vector<EffectiveNumberType> evaluate(const Vector<EffectiveNumberType>& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<EffectiveNumberType> > evaluate(const Vector< Formula<EffectiveNumberType> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<EffectiveNumberType> > evaluate(const Vector< Algebra<EffectiveNumberType> >& x) const = 0;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<EffectiveTag> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<EffectiveTag>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<EffectiveTag>* _clone() const = 0;
};

//! \relates VectorFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const VectorFunctionInterface<Void>& f) {
    return f.write(os);
}



//! \brief An interface for scalar function models on a restricted domain.
class ScalarModelInterface {
    virtual Box domain() const = 0;
    virtual ValidatedNumberType evaluate(const Vector<ValidatedNumberType>&) const = 0;
};

template<class X> class FunctionFactoryInterface;

template<> class FunctionFactoryInterface<ValidatedNumberType>
{
    typedef Box DomainType;
  public:
    virtual FunctionFactoryInterface<ValidatedTag>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunction<ValidatedTag> create(const Box& domain, const ScalarFunctionInterface<ValidatedTag>& function) const;
    inline VectorFunction<ValidatedTag> create(const Box& domain, const VectorFunctionInterface<ValidatedTag>& function) const;
    inline ScalarFunction<ValidatedTag> create_zero(const Box& domain) const;
    inline VectorFunction<ValidatedTag> create_identity(const Box& domain) const;
  private:
    virtual ScalarFunctionInterface<ValidatedTag>* _create(const Box& domain, const ScalarFunctionInterface<ValidatedTag>& function) const = 0;
    virtual VectorFunctionInterface<ValidatedTag>* _create(const Box& domain, const VectorFunctionInterface<ValidatedTag>& function) const = 0;
};

template<class X> inline OutputStream& operator<<(OutputStream& os, const FunctionFactoryInterface<ValidatedTag>& factory) {
    factory.write(os); return os;
}

} // namespace Ariadne

#endif
