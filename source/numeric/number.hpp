/***************************************************************************
 *            numeric/number.hpp
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

/*! \file numeric/number.hpp
 *  \brief Generic numbers
 */



#ifndef ARIADNE_NUMBER_HPP
#define ARIADNE_NUMBER_HPP

#include "../utility/handle.hpp"
#include "../numeric/paradigm.hpp"
#include "../utility/prototype.hpp"

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "number_interface.hpp"

#include "arithmetic.hpp"
#include "integer.hpp"
#include "dyadic.hpp"
#include "rational.hpp"
#include "real.hpp"

#include "number_interface.hpp"

namespace Ariadne {

/************ Number *********************************************************/

template<class X> struct IsNumericType;

class NumberInterface;

template<class P> class Number;
template<class P> struct IsNumericType<Number<P>> : True { };
template<class P> struct IsNumericType<LowerNumber<P>> : True { };
template<class P> struct IsNumericType<UpperNumber<P>> : True { };

struct DispatchException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

template<class PR> struct IsPrecision : False { };
template<> struct IsPrecision<DoublePrecision> : True { };
template<> struct IsPrecision<MultiplePrecision> : True { };

template<class P> Positive<Number<P>> cast_positive(Number<P> y);
template<class P> Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y);


class DeclareNumberOperators {
    friend ApproximateNumber operator+(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator-(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator*(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator/(ApproximateNumber const&, ApproximateNumber const&);

    friend ValidatedLowerNumber operator+(ValidatedLowerNumber const& y1, ValidatedLowerNumber const& y2);
    friend ValidatedLowerNumber operator-(ValidatedLowerNumber const& y1, ValidatedUpperNumber const& y2);
    friend ValidatedUpperNumber operator+(ValidatedUpperNumber const& y1, ValidatedUpperNumber const& y2);
    friend ValidatedUpperNumber operator-(ValidatedUpperNumber const& y1, ValidatedLowerNumber const& y2);

    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator+(N const& y1, D const& d2) { return y1+Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator-(N const& y1, D const& d2) { return y1-Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator*(N const& y1, D const& d2) { return y1*Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator/(N const& y1, D const& d2) { return y1/Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator+(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)+y2; }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator-(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)-y2; }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator*(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)*y2; }
    template<GenericNumber N, BuiltinFloatingPoint D> friend decltype(auto)
    operator/(D const& d1, N const& y2) { return Number<ApproximateTag>(d1)/y2; }

    template<ConcreteNumber R, class P> friend decltype(auto)
    operator+(R const& r1, Number<P> const& y2) { return Number<Paradigm<R>>(r1)+y2; }
    template<ConcreteNumber R, class P> friend decltype(auto)
    operator+(Number<P> const& y1, R const& r2) { return y1+Number<Paradigm<R>>(r2); }
};

template<class X, class P=Void> struct HasOperatorNumber {
    template<class XX, class PP, class=decltype(declval<XX>().operator Number<PP>())> static True test(int);
    template<class XX, class PP> static False test(...);
    static const bool value = decltype(test<X,P>(1))::value;
};

template<class X> struct HasOperatorNumber<X,Void> {
    template<class XX, class=decltype(declval<XX>().operator Number<Paradigm<XX>>())> static True test(int);
    template<class XX> static False test(...);
    static const bool value = decltype(test<X>(1))::value;
};


template<class X, class P> concept ConvertibleBuiltinFloatingPointToNumber
    = Same<P,ApproximateTag> and BuiltinFloatingPoint<X>;
template<class X, class P> concept ConvertibleViaRealToNumber
    = WeakerThan<P,ParadigmTag<X>> and Convertible<X,Real>;
template<class X, class P> concept ConvertibleViaNumberToNumber
    = WeakerThan<P,ParadigmTag<X>> and (not Convertible<X,Real>) and Convertible<X,Number<ParadigmTag<X>>>;

//! \ingroup NumericModule
//! \brief Generic numbers with computational paradigm \a P,  which may be %EffectiveTag, %ValidatedTag, %UpperTag, %LowerTag or %ApproximateTag.
template<class P> class Number
    : public Handle<NumberInterface>
    , public DeclareNumberOperators
{
    static_assert(IsParadigm<P>,"P must be a paradigm");
    static_assert(Same<P,ExactTag> or Same<P,EffectiveTag> or Same<P,ValidatedTag> or Same<P,ApproximateTag>);

    template<class PP> friend class Number;
    template<class PP> friend class LowerNumber;
    template<class PP> friend class UpperNumber;

    template<class PR> using ResultFloatType = FloatType<Weaker<P,ValidatedTag>,PR>;

    template<class X> static const bool IsGettableAs = IsNumericType<X>::value and IsWeaker<typename X::Paradigm,P>::value and (not Same<typename X::Paradigm,ExactTag>);
  public:
    typedef NumberInterface Interface;
    typedef P Paradigm;
    typedef Number<P> NumericType;
  public:
    using Handle<Interface>::Handle;
    explicit Number(NumberInterface* p) : Handle<NumberInterface>(p) { }
  private: public:
    Handle<NumberInterface> handle() const { return *this; }
    explicit Number(Handle<NumberInterface> h) : Handle<NumberInterface>(h) { }
  public:
    Number() : Number(Integer(0)) { }

    // Construct from a Number of a stronger paradigm
    template<StrongerThan<P> SP> Number(const Number<SP>& y) : Number<P>(y.handle()) { }

    //! Construct from a builtin integer
    template<BuiltinIntegral N> Number(const N& n) : Number<P>(Integer(n)) { }
    // Construct from a builtin floating-point number
    template<ConvertibleBuiltinFloatingPointToNumber<P> X> Number(const X& x) : Number<P>(Dyadic(x)) { }

    // Construct from a type which is convertible to Real.
    template<ConvertibleViaRealToNumber<P> X> Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    // Construct from a type which is convertible to another Number type.
    // TODO: Decide conversion properties from concrete type to Number<P>
    template<ConvertibleViaNumberToNumber<P> X>
        explicit Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    //! \brief Get the value of the number as a double-precision floating-point type
    ResultFloatType<DoublePrecision> get(DoublePrecision const& prec) const { return this->ref()._get(P(),prec); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    ResultFloatType<MultiplePrecision> get(MultiplePrecision const& prec) const { return this->ref()._get(P(),prec); }

    //! \brief Get the value of the number as a floating-point ball with the given precision and error precision.
    template<class PR, class PRE, class=EnableIf<And<IsPrecision<PR>,IsPrecision<PRE>,IsWeaker<ValidatedTag,P>>>>
    FloatBall<PR,PRE> get(PR const& prec, PRE const& errprec) const { return this->ref()._get(P(),prec,errprec); }

    //! \brief Try to dynamic_cast the object to concrete type \a X.
    template<class X> X extract() const;

    //! <p/>
    friend Number<P> operator+(Number<P> const& y) { return pos(y); } //! <p/>
    friend Number<P> operator-(Number<P> const& y) { return neg(y); } //! <p/>
    friend Number<P> operator+(Number<P> const& y1, Number<P> const& y2) { return add(y1,y2); } //! <p/>
    friend Number<P> operator-(Number<P> const& y1, Number<P> const& y2) { return sub(y1,y2); } //! <p/>
    friend Number<P> operator*(Number<P> const& y1, Number<P> const& y2) { return mul(y1,y2); } //! <p/>
    friend Number<P> operator/(Number<P> const& y1, Number<P> const& y2) { return div(y1,y2); } //! <p/>
    friend Number<P>& operator+=(Number<P>& y1, Number<P> const& y2) { return y1=y1+y2; } //! <p/>
    friend Number<P>& operator-=(Number<P>& y1, Number<P> const& y2) { return y1=y1-y2; } //! <p/>
    friend Number<P>& operator*=(Number<P>& y1, Number<P> const& y2) { return y1=y1*y2; } //! <p/>
    friend Number<P>& operator/=(Number<P>& y1, Number<P> const& y2) { return y1=y1/y2; } //! <p/>

    friend Number<P> nul(Number<P> const& y) { return _apply<P>(Nul(),y); } //! <p/>
    friend Number<P> pos(Number<P> const& y) { return _apply<P>(Pos(),y); } //! <p/>
    friend Number<P> neg(Number<P> const& y) { return _apply<P>(Neg(),y); } //! <p/>
    friend Number<P> sqr(Number<P> const& y) { return _apply<P>(Sqr(),y); } //! <p/>
    friend Number<P> hlf(Number<P> const& y) { return _apply<P>(Hlf(),y); } //! <p/>
    friend Number<P> rec(Number<P> const& y) { return _apply<P>(Rec(),y); } //! <p/>
    friend Number<P> add(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Add(),y1,y2); } //! <p/>
    friend Number<P> sub(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Sub(),y1,y2); } //! <p/>
    friend Number<P> mul(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Mul(),y1,y2); } //! <p/>
    friend Number<P> div(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Div(),y1,y2); } //! <p/>
    friend Number<P> sqrt(Number<P> const& y) { return _apply<P>(Sqrt(),y); } //! <p/>
    friend Number<P> exp(Number<P> const& y) { return _apply<P>(Exp(),y); } //! <p/>
    friend Number<P> log(Number<P> const& y) { return _apply<P>(Log(),y); } //! <p/>
    friend Number<P> sin(Number<P> const& y) { return _apply<P>(Sin(),y); } //! <p/>
    friend Number<P> cos(Number<P> const& y) { return _apply<P>(Cos(),y); } //! <p/>
    friend Number<P> tan(Number<P> const& y) { return _apply<P>(Tan(),y); } //! <p/>
    friend Number<P> asin(Number<P> const& y) { return _apply<P>(Asin(),y); } //! <p/>
    friend Number<P> acos(Number<P> const& y) { return _apply<P>(Acos(),y); } //! <p/>
    friend Number<P> atan(Number<P> const& y) { return _apply<P>(Atan(),y); } //! <p/>

    friend Number<P> pow(Number<P> const& y, Nat m) { return _apply<P>(Pow(),y,m); } //! <p/>
    friend Number<P> pow(Number<P> const& y, Int n) { return _apply<P>(Pow(),y,n); } //! <p/>

    friend Positive<Number<P>> abs(Number<P> const& y) { return Positive<Number<P>>(_apply<P>(Abs(),y)); } //! <p/>
    friend Number<P> max(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Max(),y1,y2); } //! <p/>
    friend Number<P> min(Number<P> const& y1, Number<P> const& y2) { return _apply<P>(Min(),y1,y2); } //! <p/>


    friend LogicalType<Equality<P>> operator==(Number<P> const& y1, Number<P> const& y2) {
        return LogicalType<Equality<P>>(y1.ref()._equals(y2.ref())); } //! <p/>
    friend LogicalType<LessThan<P>> operator< (Number<P> const& y1, Number<P> const& y2) {
        return LogicalType<LessThan<P>>(y1.ref()._less(y2.ref())); } //! <p/>
    friend LogicalType<LessThan<P>> operator> (Number<P> const& y1, Number<P> const& y2) { return (y2<y1); } //! <p/>
    friend LogicalNegationType<LogicalType<Equality<P>>> operator!=(Number<P> const& y1, Number<P> const& y2) { return !(y1==y2); } //! <p/>
    friend LogicalType<LessThan<P>> operator<=(Number<P> const& y1, Number<P> const& y2) { return !(y1>y2); } //! <p/>
    friend LogicalType<LessThan<P>> operator>=(Number<P> const& y1, Number<P> const& y2) { return !(y1<y2); }

    //! \brief The name of the dynamic type held by the %Number object.
    String class_name() const { return this->ref()._class_name(); }

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, Number<P> const& y) { return y.ref()._write(os); }\
  private:
    template<class RP, class OP> static Number<RP> _apply(OP op, Number<P> const& y) { return Number<RP>(y.ref()._apply(op)); }
    template<class RP, class OP> static Number<RP> _apply(OP op, Number<P> const& y, Int n) { return Number<RP>(y.ref()._apply(op,n)); }
    template<class RP, class OP, class P1, class P2> static Number<RP> _apply(OP op, Number<P1> const& y1, Number<P2> const& y2) {
        return Number<RP>(y1.ref()._apply(op,&y2.ref())); }
};


//! \ingroup NumericModule
//! \brief Generic lower (real) numbers with computational paradigm \a P, which may be %EffectiveTag or %ValidatedTag.
template<class P> class LowerNumber
{
    static_assert(IsSame<P,EffectiveTag>::value or IsSame<P,ValidatedTag>::value,"P must be a paradigm");
    friend class UpperNumber<P>;
  private: public:
    Handle<NumberInterface> _handle;
    explicit LowerNumber(NumberInterface* p) : _handle(p) { }
  private: public:
    explicit LowerNumber(Handle<NumberInterface> h) : _handle(h) { }
    Handle<NumberInterface> handle() const { return this->_handle; }
  private:
    NumberInterface const& ref() const { return this->_handle.reference(); }
  public:
    typedef P Paradigm;
    typedef LowerNumber<P> NumericType;

    LowerNumber() : LowerNumber(Integer(0)) { }

    //! \brief Construct from a LowerNumber of a stronger paradigm
    template<StrongerThan<P> SP> LowerNumber(const LowerNumber<SP>& y) : LowerNumber<P>(y.handle()) { }
    //! \brief Construct from a type convertible to a Number.
    template<ConvertibleTo<Number<P>> X> LowerNumber(const X& x) : LowerNumber<P>(Number<P>(x).handle()) { }

    template<class PR> FloatLowerBound<PR> get(PR pr) const { return this->ref()._get(LowerTag(),pr); }

    template<class X> X extract() const;

    friend LowerNumber<P> operator+(LowerNumber<P> const& y) { return pos(y); }
    friend UpperNumber<P> operator-(LowerNumber<P> const& y) { return neg(y); }
    friend LowerNumber<P> operator+(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return add(y1,y2); }
    friend LowerNumber<P> operator-(LowerNumber<P> const& y1, UpperNumber<P> const& y2) { return sub(y1,y2); }
    friend LowerNumber<P>& operator+=(LowerNumber<P>& y1, LowerNumber<P> const& y2) { return y1=y1+y2; }
    friend LowerNumber<P>& operator-=(LowerNumber<P>& y1, UpperNumber<P> const& y2) { return y1=y1-y2; }

    friend LowerNumber<P> pos(LowerNumber<P> const& y) { return LowerNumber<P>(y.ref()._apply(Pos())); }
    friend UpperNumber<P> neg(LowerNumber<P> const& y) { return UpperNumber<P>(y.ref()._apply(Neg())); }
    friend LowerNumber<P> add(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return LowerNumber<P>(y1.ref()._apply(Add(),&y2.ref())); }
    friend LowerNumber<P> sub(LowerNumber<P> const& y1, UpperNumber<P> const& y2) { return LowerNumber<P>(y1.ref()._apply(Sub(),&y2.handle().reference())); }

    friend LowerNumber<P> sqrt(LowerNumber<P> const& y) { return LowerNumber<P>(y.ref()._apply(Sqrt())); }
    friend LowerNumber<P> exp(LowerNumber<P> const& y) { return LowerNumber<P>(y.ref()._apply(Exp())); }
    friend LowerNumber<P> log(LowerNumber<P> const& y) { return LowerNumber<P>(y.ref()._apply(Log())); }
    friend LowerNumber<P> atan(LowerNumber<P> const& y) { return LowerNumber<P>(y.ref()._apply(Atan())); }

    friend LowerNumber<P> max(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return LowerNumber<P>(y1.ref()._apply(Max(),&y2.ref())); }
    friend LowerNumber<P> min(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return LowerNumber<P>(y1.ref()._apply(Min(),&y2.ref())); }

    friend UpperLogicalType<P> operator==(LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return UpperLogicalType<P>(y1.ref()._equals(y2.handle().reference())); }
    friend LowerLogicalType<P> operator!=(LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return not (y1 == y2); }
    friend UpperLogicalType<P> operator< (LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return UpperLogicalType<P>(y1.ref()._less(y2.handle().reference())); }
    friend LowerLogicalType<P> operator> (LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return y2 <  y1; }

    String class_name() const { return this->ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, LowerNumber<P> const& y) { return y.ref()._write(os); }
};


//! \ingroup NumericModule
//! \brief Generic upper (real) numbers with computational paradigm \a P, which may be %EffectiveTag or %ValidatedTag.
template<class P> class UpperNumber
{
    static_assert(IsSame<P,EffectiveTag>::value or IsSame<P,ValidatedTag>::value,"P must be a paradigm");
    friend class LowerNumber<P>;
  private: public:
    Handle<NumberInterface> _handle;
    explicit UpperNumber(NumberInterface* p) : _handle(p) { }
  private: public:
    explicit UpperNumber(Handle<NumberInterface> h) : _handle(h) { }
    Handle<NumberInterface> handle() const { return this->_handle; }
  private:
    NumberInterface const& ref() const { return this->_handle.reference(); }
  public:
    typedef P Paradigm;
    typedef UpperNumber<P> NumericType;

    UpperNumber() : UpperNumber(Integer(0)) { }

    //! \brief Construct from a UpperNumber of a stronger paradigm
    template<StrongerThan<P> SP> UpperNumber(const UpperNumber<SP>& y) : UpperNumber<P>(y.handle()) { }
    //! \brief Construct from a type convertible to a Number.
    template<ConvertibleTo<Number<P>> X> UpperNumber(const X& x) : UpperNumber<P>(Number<P>(x).handle()) { }

    template<class PR> FloatUpperBound<PR> get(PR pr) const { return this->ref()._get(UpperTag(),pr); }

    template<class X> X extract() const;

    friend UpperNumber<P> operator+(UpperNumber<P> const& y) { return pos(y); }
    friend LowerNumber<P> operator-(UpperNumber<P> const& y) { return neg(y); }
    friend UpperNumber<P> operator+(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return add(y1,y2); }
    friend UpperNumber<P> operator-(UpperNumber<P> const& y1, LowerNumber<P> const& y2) { return sub(y1,y2); }
    friend UpperNumber<P>& operator+=(UpperNumber<P>& y1, UpperNumber<P> const& y2) { return y1=y1+y2; }
    friend UpperNumber<P>& operator-=(UpperNumber<P>& y1, LowerNumber<P> const& y2) { return y1=y1-y2; }

    friend UpperNumber<P> pos(UpperNumber<P> const& y) { return UpperNumber<P>(y.ref()._apply(Pos())); }
    friend LowerNumber<P> neg(UpperNumber<P> const& y) { return LowerNumber<P>(y.ref()._apply(Neg())); }
    friend UpperNumber<P> add(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return UpperNumber<P>(y1.ref()._apply(Add(),&y2.ref())); }
    friend UpperNumber<P> sub(UpperNumber<P> const& y1, LowerNumber<P> const& y2) { return UpperNumber<P>(y1.ref()._apply(Sub(),&y2.handle().reference())); }

    friend UpperNumber<P> sqrt(UpperNumber<P> const& y) { return UpperNumber<P>(y.ref()._apply(Sqrt())); }
    friend UpperNumber<P> exp(UpperNumber<P> const& y) { return UpperNumber<P>(y.ref()._apply(Exp())); }
    friend UpperNumber<P> log(UpperNumber<P> const& y) { return UpperNumber<P>(y.ref()._apply(Log())); }
    friend UpperNumber<P> atan(UpperNumber<P> const& y) { return UpperNumber<P>(y.ref()._apply(Atan())); }

    friend UpperNumber<P> max(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return UpperNumber<P>(y1.ref()._apply(Max(),&y2.ref())); }
    friend UpperNumber<P> min(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return UpperNumber<P>(y1.ref()._apply(Min(),&y2.ref())); }

    friend UpperLogicalType<P> operator==(UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        return UpperLogicalType<P>(y1.ref()._equals(y2.handle().reference())); }
    friend LowerLogicalType<P> operator!=(UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        return not (y1 == y2); }
    friend LowerLogicalType<P> operator< (UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        return LowerLogicalType<P>(y1.ref()._less(y2.handle().reference())); }
    friend UpperLogicalType<P> operator> (UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        return y2 <  y1; }

    String class_name() const { return this->ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, UpperNumber<P> const& y) { return y.ref()._write(os); }
};



template<class P> class Positive<Number<P>> : public Number<P> {
    friend Number<P> const& unsign(Positive<Number<P>> const& y) { return y; }
  public:
    Positive<Number<P>>() : Number<P>() { }
    explicit Positive<Number<P>>(Number<P> const& y) : Number<P>(y) { }
    template<BuiltinUnsignedIntegral N> Positive<Number<P>>(N n) : Number<P>(n) { }
    template<class N> requires Convertible<N,ExactNumber>
        Positive<Number<P>>(const Positive<N>& n) : Number<P>(ExactNumber(static_cast<N const&>(n))) { }
    template<class N> requires Constructible<ExactNumber,N> and (not BuiltinIntegral<N>)
        explicit Positive<Number<P>>(const N& n) : Number<P>(ExactNumber(n)) { }
    explicit operator Number<P> () const { return *this; }

    friend Positive<Number<P>> operator+(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> operator*(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) {
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> add(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> mul(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) {
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> max(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) {
        return cast_positive(max(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> abs(Positive<Number<P>> const& y) { return y; }
};

template<class P> class Positive<LowerNumber<P>> : public LowerNumber<P> {
    friend LowerNumber<P> const& unsign(Positive<LowerNumber<P>> const& y) { return y; }
  public:
    Positive<LowerNumber<P>>() : LowerNumber<P>() { }
    explicit Positive<LowerNumber<P>>(LowerNumber<P> const& y) : LowerNumber<P>(y) { }
    template<BuiltinUnsignedIntegral N>
        Positive<LowerNumber<P>>(N n) : LowerNumber<P>(n) { }
    template<class N> requires Constructible<ExactNumber,N>
        Positive<LowerNumber<P>>(const Positive<N>& n) : LowerNumber<P>(ExactNumber(static_cast<N const&>(n))) { }
    template<class N> requires Constructible<ExactNumber,N> and (not BuiltinIntegral<N>)
        Positive<LowerNumber<P>>(const N& n) : LowerNumber<P>(ExactNumber(n)) { }
    explicit operator LowerNumber<P> () const { return *this; }

    friend LowerNumber<P> mul(LowerNumber<P> const& y1, LowerNumber<P> const& y2);

    friend Positive<LowerNumber<P>> operator+(Positive<LowerNumber<P>> const& y1, Positive<LowerNumber<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<LowerNumber<P>> operator*(Positive<LowerNumber<P>> const& y1, Positive<LowerNumber<P>> const& y2) {
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<LowerNumber<P>> add(Positive<LowerNumber<P>> const& y1, Positive<LowerNumber<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<LowerNumber<P>> mul(Positive<LowerNumber<P>> const& y1, Positive<LowerNumber<P>> const& y2) {
        return cast_positive(Number<P>(y1.ref()._apply(Mul(),y2.ref()))); }
    friend Positive<LowerNumber<P>> max(Positive<LowerNumber<P>> const& y1, Positive<LowerNumber<P>> const& y2) {
        return cast_positive(max(unsign(y1),unsign(y2))); }
};

template<class P> class Positive<UpperNumber<P>> : public UpperNumber<P> {
    friend UpperNumber<P> const& unsign(Positive<UpperNumber<P>> const& y) { return y; }
  public:
    Positive<UpperNumber<P>>() : UpperNumber<P>() { }
    explicit Positive<UpperNumber<P>>(UpperNumber<P> const& y) : UpperNumber<P>(y) { }
    template<BuiltinUnsignedIntegral N>
        Positive<UpperNumber<P>>(N n) : UpperNumber<P>(n) { }
    template<class N> requires Constructible<ExactNumber,N>
        Positive<UpperNumber<P>>(const Positive<N>& n) : UpperNumber<P>(ExactNumber(static_cast<N const&>(n))) { }
    template<class N> requires Constructible<ExactNumber,N> and (not BuiltinIntegral<N>)
        Positive<UpperNumber<P>>(const N& n) : UpperNumber<P>(ExactNumber(n)) { }
    explicit operator UpperNumber<P> () const { return *this; }

    friend UpperNumber<P> mul(UpperNumber<P> const& y1, UpperNumber<P> const& y2);

    friend Positive<UpperNumber<P>> operator+(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<UpperNumber<P>> operator*(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<UpperNumber<P>> add(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<UpperNumber<P>> mul(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(Number<P>(y1.ref()._apply(Mul(),y2.ref()))); }
    friend Positive<UpperNumber<P>> max(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(max(unsign(y1),unsign(y2))); }
};

template<class X> decltype(auto) cast_generic(X const& x) { return x.generic(); }

template<class P> Positive<Number<P>> cast_positive(Number<P> y) { return Positive<Number<P>>(y); }
template<class P> Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y) { return Positive<UpperNumber<P>>(y); }

template<class P> Positive<ExactNumber> cast_exact(Positive<Number<P>> const& y) { return cast_positive(cast_exact(unsign(y))); }
template<class P> Positive<ExactNumber> cast_exact(Positive<UpperNumber<P>> const& y) { return cast_positive(cast_exact(unsign(y))); }


} // namespace Ariadne

#endif
