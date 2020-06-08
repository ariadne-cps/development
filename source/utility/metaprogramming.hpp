/***************************************************************************
 *            utility/metaprogramming.hpp
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

/*! \file utility/metaprogramming.hpp
 *  \brief Classes for template metaprogramming.
 */

#ifndef ARIADNE_METAPROGRAMMING_HPP
#define ARIADNE_METAPROGRAMMING_HPP

#include <type_traits>
#include <concepts>

namespace Ariadne {

typedef void Void;
typedef bool Bool;
typedef std::size_t SizeType;

using std::declval;

using True = std::true_type;
using False = std::false_type;

template<class T, T V> using IntegralConstant = std::integral_constant<T,V>;

template<bool b, class... PS> struct ValueAnd;
template<class... PS> struct And;
template<class P1, class... PS> struct And<P1,PS...> : ValueAnd<P1::value,PS...> { };
template<> struct And<> : True { };
template<class... PS> struct ValueAnd<false, PS...> : False { };
template<class... PS> struct ValueAnd<true, PS...> : And<PS...> { };

template<class P> struct Not { static const bool value = !P::value; };

template<class T1, class... TS> struct First { typedef T1 Type; };
template<class T1, class T2, class... TS> struct Second { typedef T2 Type; };

typedef bool Dummy;

template<class P, typename V=Dummy> using EnableIf = typename std::enable_if<P::value,V>::type;
template<class P, typename V=Dummy> using DisableIf = typename std::enable_if<not P::value,V>::type;


template<class T1, class T2> using IsSame = std::is_same<T1,T2>;
template<class T, class U> using IsConvertible = std::is_convertible<T,U>;
template<class T> using IsBuiltinIntegral = std::is_integral<T>;
template<class T> using IsBuiltinFloatingPoint = std::is_floating_point<T>;
/*
template<class T, T V> using IntegralConstant = std::integral_constant<T,V>;
template<class T1, class T2> using IsSame = std::is_same<T1,T2>;
template<class T, class U> using IsConvertible = std::is_convertible<T,U>;
template<class T, class... U> using IsConstructible = std::is_constructible<T,U...>;
template<class T, class U> using IsStrictlyConstructible = And<IsConstructible<T,U>,Not<IsConvertible<U,T>>>;
template<class T> using IsDefaultConstructible = std::is_default_constructible<T>;
template<class T> using IsCopyConstructible = std::is_copy_constructible<T>;
template<class T, class U> using IsAssignable = std::is_assignable<T,U>;
template<class T, class U> using IsBaseOf = std::is_base_of<T,U>;
//template<class F, class... AS> using IsInvocable = std::is_invokable<F,AS...>;
//template<class R, class F, class... AS> using IsInvocableReturning = std::is_invokable_r<R,F,AS...>;
template<class T> using IsBuiltinArithmetic = std::is_arithmetic<T>;
template<class T> using IsBuiltinIntegral = std::is_integral<T>;
template<class T> using IsBuiltinFloatingPoint = std::is_floating_point<T>;
template<class T> using IsBuiltinSigned = std::is_signed<T>;
template<class T> using IsBuiltinUnsigned = std::is_unsigned<T>;
template<class T> using IsBuiltinSignedIntegral = std::integral_constant<bool,std::is_integral<T>::value and std::is_signed<T>::value>;
template<class T> using IsBuiltinUnsignedIntegral = std::integral_constant<bool,std::is_integral<T>::value and std::is_unsigned<T>::value>;
*/

template<class T, class... US> struct IsOneOf;
template<class T> struct IsOneOf<T> : False { };
template<class T, class... US> struct IsOneOf<T,T,US...> : True { };
template<class T, class U0, class... US> struct IsOneOf<T,U0,US...> : IsOneOf<T,US...> { };

template<class T, class... TS> struct IndexOf;
template<class T, class... TS> struct IndexOf<T,T,TS...> { static const SizeType N=0; };
template<class T, class T0, class... TS> struct IndexOf<T,T0,TS...> { static const SizeType N=IndexOf<T,TS...>::N+1u; };
template<class T> struct IndexOf<T> { };


template<class F, class... AS> using InvokeResult = typename std::invoke_result<F,AS...>::type;
template<class SIG> struct ResultOfTrait;
template<class F, class... AS> struct ResultOfTrait<F(AS...)> { typedef typename std::invoke_result<F,AS...>::type Type; };
template<class SIG> using ResultOf = typename ResultOfTrait<SIG>::Type;

template<class T> struct Self { typedef T Type; };
template<class T> using SelfType = typename Self<T>::Type;


struct Fallback { };
struct DontCare { template<class T> DontCare(T); };

template<class R> struct Return { typedef R ReturnType; };
template<class RET> using ReturnType = typename RET::ReturnType;

template<class F, class A, class R=DontCare, class=Fallback> struct IsCallable : False { };
template<class F, class A, class R> struct IsCallable<F,A,R, EnableIf<IsConvertible<decltype(declval<F>()(declval<A>())),R>,Fallback>> : True { };


template<class X> using RemoveConst = typename std::remove_const<X>::type;
template<class X> using RemoveReference = typename std::remove_reference<X>::type;

template<class L> using LogicalNegationType = decltype(!declval<L>());
template<class L1, class L2=L1> using LogicalConjunctionType = decltype(declval<L1>() && declval<L2>());
template<class L1, class L2=L1> using LogicalDisjunctionType = decltype(declval<L1>() || declval<L2>());

template<class X> using NegationType = decltype(-declval<X>());
template<class X1, class X2=X1> using SumType = decltype(declval<X1>()+declval<X2>());
template<class X1, class X2=X1> using DifferenceType = decltype(declval<X1>()-declval<X2>());
template<class X1, class X2=X1> using ProductType = decltype(declval<X1>()*declval<X2>());
template<class X1, class X2=X1> using QuotientType = decltype(declval<X1>()/declval<X2>());
template<class X1, class X2=X1> using ArithmeticType = SumType<ProductType<X1,X2>>;
template<class X1, class X2=X1> using EqualsType = decltype(declval<X1>()==declval<X2>());
template<class X1, class X2=X1> using LessType = decltype(declval<X1>()< declval<X2>());
template<class X1, class X2=X1> using ComparisonType = LessType<X1,X2>;

template<class X> using MagType = decltype(mag(declval<X>()));

template<class X1, class X2=X1> using EqualityType = decltype(declval<X1>()==declval<X2>());
template<class X1, class X2=X1> using InequalityType = decltype(declval<X1>()!=declval<X2>());

template<class X1, class X2=X1> using InplaceSumType = RemoveReference<decltype(declval<X1&>()+=declval<X2>())>;
template<class X1, class X2=X1> using InplaceDifferenceType = RemoveReference<decltype(declval<X1&>()-=declval<X2>())>;
template<class X1, class X2=X1> using InplaceProductType = RemoveReference<decltype(declval<X1&>()*=declval<X2>())>;
template<class X1, class X2=X1> using InplaceQuotientType = RemoveReference<decltype(declval<X1&>()/=declval<X2>())>;


template<class SIG> using ResultOf = typename std::result_of<SIG>::type;
template<class F, class... AS> using InvokeResult = typename std::invoke_result<F,AS...>::type;

template<class F, class... AS> struct IsInvocable;
template<class R, class F, class... AS> struct IsInvocableReturning;

template<class F, class A> struct IsInvocable<F,A> {
    template<class FF, class AA, class=decltype(std::declval<FF>()(std::declval<AA>()))>
        static std::true_type test(int);
    template<class FF, class AA>
        static std::false_type test(...);
    static const bool value = decltype(test<F,A>(1))::value;
};

template<class R, class F, class A> struct IsInvocableReturning<R,F,A> {
    template<class RR, class FF, class AA, class=decltype(std::declval<RR>()=std::declval<FF>()(std::declval<AA>()))>
        static std::true_type test(int);
    template<class RR, class FF, class AA>
        static std::false_type test(...);
    static const bool value = decltype(test<R,F,A>(1))::value;
};




template<class T1, class T2> concept SameAs = std::same_as<T1,T2>;
template<class T1, class T2> concept Same = std::same_as<T1,T2>;
template<class T> concept BuiltinArithmetic = std::is_arithmetic<T>::value;
template<class T> concept BuiltinIntegral = std::integral<T>;
template<class T> concept BuiltinFloatingPoint = std::floating_point<T>;
template<class T> concept BuiltinUnsignedIntegral = std::unsigned_integral<T>;
template<class T> concept BuiltinSignedIntegral = std::signed_integral<T>;
//template<class F, class T> concept ConvertibleTo = std::convertible_to<F,T>;
//template<class F, class T> concept Convertible = std::convertible_to<F,T>;
template<class F, class T> concept Convertible = std::is_convertible<F,T>::value;
template<class F, class T> concept ConvertibleTo = std::is_convertible<F,T>::value;
template<class T, class... US> concept ConstructibleFrom = std::constructible_from<T,US...>;
template<class T, class... US> concept Constructible = std::constructible_from<T,US...>;
template<class T> concept DefaultInitializable = std::default_initializable<T>;
template<class T> concept DefaultConstructible = std::default_initializable<T>;
template<class T> concept CopyConstructible = std::copy_constructible<T>;
//template<class F, class T> concept ExplicitlyConvertibleTo = std::constructible_from<T,F> and not std::convertible_to<F,T>;
//template<class F, class T> concept ExplicitlyConvertible = std::constructible_from<T,F> and not std::convertible_to<F,T>;
template<class F, class T> concept ExplicitlyConvertibleTo = std::is_constructible<T,F>::value and not std::is_convertible<F,T>::value;
template<class F, class T> concept ExplicitlyConvertible = std::is_constructible<T,F>::value and not std::is_convertible<F,T>::value;

//template<class T, class F> concept AssignableFrom = std::assignable_from<T,F>;
//template<class F, class T> concept AssignableTo = std::assignable_from<T,F>;
//template<class T, class F> concept Assignable = std::assignable_from<T,F>;
template<class T, class F> concept Assignable = std::is_assignable<T,F>::value;
template<class T, class F> concept AssignableFrom = std::is_assignable<T,F>::value;
template<class F, class T> concept AssignableTo = std::is_assignable<T,F>::value;
template<class D, class B> concept DerivedFrom = std::derived_from<D,B>;
template<class B, class D> concept BaseOf = std::derived_from<D,B>;
template<class F, class... AS> concept Invocable = std::invocable<F,AS...>;

template<class R, class F, class... AS> concept InvocableReturning = IsInvocableReturning<R,F,AS...>::value;
template<class T, class... TS> concept OneOf = IsOneOf<T,TS...>::value;
template<class T1, class T2, class T3> concept AreSame = Same<T1,T2> and Same<T2,T3>;


template<class A1, class A2> concept HasEquality = requires(A1 a1, A2 a2) { { a1==a2 }; };
template<class A1, class A2> concept CanAdd = requires(A1 a1, A2 a2) { { a1+a2 }; };
template<class A1, class A2> concept CanSubtract = requires(A1 a1, A2 a2) { { a1-a2 }; };
template<class A1, class A2> concept CanMultiply = requires(A1 a1, A2 a2) { { a1*a2 }; };
template<class A1, class A2> concept CanDivide = requires(A1 a1, A2 a2) { { a1/a2 }; };


} // namespace Ariadne

#endif
