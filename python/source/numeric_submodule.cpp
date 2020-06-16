/***************************************************************************
 *            numeric_submodule.cpp
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

#include "pybind11.hpp"
#include <pybind11/operators.h>

#include "utilities.hpp"

#if defined(__GNUG__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wattributes"
#endif

#include "numeric/logical.hpp"
#include "numeric/accuracy.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/float.hpp"
#include "numeric/complex.hpp"

#include "numeric/casts.hpp"


namespace Ariadne {

static const bool ALLOW_CONCRETE_TO_GENERIC_NUMBER_CONVERSIONS = true;
//static const bool ALLOW_CONCRETE_TO_GENERIC_NUMBER_CONVERSIONS = false;
//static const bool ALLOW_CONCRETE_TO_GENERIC_NUMBER_CONVERSIONS = IsConvertible<FloatDPBounds,ValidatedNumber>::value;

extern template Nat Error<FloatDP>::output_places;
extern template Nat Error<FloatMP>::output_places;

template<> String class_name<DoublePrecision>() { return "DoublePrecision"; }
template<> String class_name<MultiplePrecision>() { return "MultiplePrecision"; }

template<class T> String numeric_class_tag();
template<> String numeric_class_tag<DoublePrecision>() { return "DP"; }
template<> String numeric_class_tag<MultiplePrecision>() { return "MP"; }

template<class T> struct PythonName { const char* get() const { return class_name<T>().c_str(); } };
template<class T> inline const char* python_name() { return PythonName<T>().get(); }

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Integer>& repr) {
    return os << "Integer("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Dyadic>& repr) {
    return os << "Dyadic("<<repr.reference().mantissa()<<","<<repr.reference().exponent()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr) {
    return os << "Rational("<<repr.reference().numerator()<<","<<repr.reference().denominator()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr) {
    return os << "Real("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<DyadicBounds>& repr) {
    return os << "DyadicBounds("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RationalBounds>& repr) {
    return os << "RationalBounds("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloatDP>& repr) {
    return os << "FloatDP("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloatMP>& repr) {
    return os << "FloatMP("<<repr.reference()<<")"; }

template<class P> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Number<P>>& repr) {
    return os << class_name<Number<P>>()<<"("<<repr.reference()<<")"; }
template<class P> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperNumber<P>>& repr) {
    return os << class_name<UpperNumber<P>>()<<"("<<repr.reference()<<")"; }
template<class P> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<LowerNumber<P>>& repr) {
    return os << class_name<LowerNumber<P>>()<<"("<<repr.reference()<<")"; }

template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Approximation<F>>& repr) {
    return os << class_name<F>() << "Approximation("<<repr.reference().raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<LowerBound<F>>& repr) {
    return os << class_name<F>() << "LowerBound("<<repr.reference().raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperBound<F>>& repr) {
    return os << class_name<F>() << "UpperBound("<<repr.reference().raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Bounds<F>>& repr) {
    return os << class_name<F>() << "Bounds("<<repr.reference().lower().raw()<<","<<repr.reference().upper().raw()<<")"; }
template<class F, class FE> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Ball<F,FE>>& repr) {
    return os << class_name<F>() << "Ball("<<repr.reference().value_raw()<<","<<repr.reference().error_raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Value<F>>& repr) {
    return os << class_name<F>() << "Value("<<repr.reference().raw()<<")"; }
template<class FE> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Error<FE>>& repr) {
    return os << class_name<FE>() << "Error("<<repr.reference().raw()<<")"; }



template<class L> Bool _decide_(L l) { return decide(l); }
template<class L> Bool _definitely_(L l) { return definitely(l); }
template<class L> Bool _possibly_(L l) { return possibly(l); }
template<class L> Bool _probably_(L l) { return probably(l); }

template<class L, class E=Effort> auto _check_(L const& l, E const& e) -> decltype(l.check(e)) { return l.check(e); }

template<class OP, class... TS> auto py_apply(TS const& ... ts) -> decltype(OP()(ts...)){ OP op; return op(ts...); }

template<class A> auto _is_nan_(A const& a) -> Bool { return is_nan(a); }
template<class A> auto _is_inf_(A const& a) -> Bool { return is_inf(a); }
template<class A> auto _is_finite_(A const& a) -> Bool { return is_finite(a); }
template<class A> auto _is_zero_(A const& a) -> Bool { return is_zero(a); }

template<class X> Void define_infinitary_checks(pybind11::module& module, pybind11::class_<X>& pyclass) {
    module.def("is_nan", &_is_nan_<X>);
    module.def("is_inf", &_is_inf_<X>);
    module.def("is_finite", &_is_finite_<X>);
    module.def("is_zero", &_is_zero_<X>);
}

template<class X> Void define_infinitary(pybind11::module& module, pybind11::class_<X>& pyclass) {
    if constexpr (HasPrecisionType<X>) {
        typedef typename X::PrecisionType PR;
        pyclass.def_static("nan", (X(*)(PR)) &X::nan);
        pyclass.def_static("inf", (X(*)(PR)) &X::inf);
        pyclass.def_static("inf", (X(*)(Sign,PR)) &X::inf);
    } else {
        pyclass.def_static("nan", (X(*)()) &X::nan);
        pyclass.def_static("inf", (X(*)()) &X::inf);
        pyclass.def_static("inf", (X(*)(Sign)) &X::inf);
    }
    define_infinitary_checks(module,pyclass);
}

template<class L, class U> Pair<L,U> pair_from_dict(pybind11::dict dct) {
    assert(dct.size()==1);
    pybind11::detail::dict_iterator::reference item = *dct.begin();
    pybind11::handle lh = item.first;
    pybind11::handle uh = item.second;
    return Pair<L,U>(pybind11::cast<L>(lh),pybind11::cast<U>(uh));
}

template<class X> X bounds_from_dict(pybind11::dict dct) {
    typedef decltype(declval<X>().lower_raw()) L;
    typedef decltype(declval<X>().upper_raw()) U;
    assert(dct.size()==1);
    pybind11::detail::dict_iterator::reference item = *dct.begin();
    pybind11::handle lh = item.first;
    pybind11::handle uh = item.second;
    L l = pybind11::cast<L>(lh);
    U u = pybind11::cast<U>(uh);
    return X(l,u);
}

template<class X, class PR> X bounds_from_dict_and_properties(pybind11::dict dct, PR pr) {
    typedef ValidatedLowerNumber L;
    typedef ValidatedUpperNumber U;
    assert(dct.size()==1);
    pybind11::detail::dict_iterator::reference item = *dct.begin();
    pybind11::handle lh = item.first;
    pybind11::handle uh = item.second;
    L l = pybind11::cast<L>(lh);
    U u = pybind11::cast<U>(uh);
    return X(l,u,pr);
}


} // namespace Ariadne


using namespace Ariadne;
using pymodule = pybind11::module;
using pybind11::init;
using pybind11::detail::self;
using pybind11::implicitly_convertible;


Void export_effort(pymodule& module) {
    pybind11::class_<Effort> effort_class(module,"Effort");
    effort_class.def(init<Nat>());
    effort_class.def("work",&Effort::work);
    effort_class.def("__str__", &__cstr__<Effort>);
}

template<class L> void export_effective_logical(pymodule& module, std::string name)
{
    typedef decltype(declval<L>().check(declval<Effort>())) CheckType;
    OutputStream& operator<<(OutputStream& os, L l);

    pybind11::class_<L> logical_class(module,name.c_str());
    logical_class.def(init<bool>());
    logical_class.def(init<L>());
    logical_class.def("__str__", &__cstr__<L>);
    logical_class.def("__repr__", &__cstr__<L>);
    logical_class.def("check", (CheckType(L::*)(Effort)) &L::check);
    module.def("check", [](L l, Effort eff){return check(l,eff);});
    define_logical(module,logical_class);
}

template<class L> void export_logical(pymodule& module, std::string name)
{
    OutputStream& operator<<(OutputStream& os, L l);
    pybind11::class_<L> logical_class(module,name.c_str());
    logical_class.def(init<bool>());
    logical_class.def(init<L>());
    logical_class.def("__str__", &__cstr__<L>);
    logical_class.def("__repr__", &__cstr__<L>);
    define_logical(module,logical_class);

    module.def("decide", &_decide_<L>);
    module.def("possibly", &_possibly_<L>);
    module.def("definitely", &_definitely_<L>);

}

template<> void export_logical<Boolean>(pymodule& module, std::string name)
{
    typedef Boolean L;
    OutputStream& operator<<(OutputStream& os, L l);
    pybind11::class_<L> logical_class(module,name.c_str());
    logical_class.def(init<bool>());
    logical_class.def(init<L>());
    logical_class.def("__str__", &__cstr__<L>);
    logical_class.def("__repr__", &__cstr__<L>);
    define_logical(module,logical_class);

//    implicitly_convertible<LogicalType<ExactTag>,bool>();
}



Void export_logicals(pymodule& module) {
    export_logical<Boolean>(module,"Boolean");
    export_effective_logical<Sierpinskian>(module,"Sierpinskian");
    export_effective_logical<NegatedSierpinskian>(module,"NegatedSierpinskian");
    export_effective_logical<Kleenean>(module,"Kleenean");
    export_effective_logical<LowerKleenean>(module,"LowerKleenean");
    export_effective_logical<UpperKleenean>(module,"UpperKleenean");
    export_logical<ValidatedKleenean>(module,"ValidatedKleenean");
    export_logical<ValidatedUpperKleenean>(module,"ValidatedUpperKleenean");
    export_logical<ValidatedLowerKleenean>(module,"ValidatedLowerKleenean");
    export_logical<ApproximateKleenean>(module,"ApproximateKleenean");
}


Void export_accuracy(pymodule& module) {
    pybind11::class_<Accuracy> accuracy_class(module,"Accuracy");
    accuracy_class.def(init<Dyadic>());
    accuracy_class.def("error",&Accuracy::error);
    accuracy_class.def("__str__", &__cstr__<Accuracy>);
    accuracy_class.def("__repr__", &__cstr__<Accuracy>);
}



void export_builtins(pymodule& module)
{
    pybind11::class_<ApproximateDouble> approximate_double_class(module,"ApproximateDouble");
    approximate_double_class.def(init<double>());
    approximate_double_class.def("__neg__", &__neg__<ApproximateDouble>);
    approximate_double_class.def("__str__", &__cstr__<ApproximateDouble>);
    implicitly_convertible<double,ApproximateDouble>();
    approximate_double_class.def("__hash__", [](ApproximateDouble const& ad){std::hash<double> hasher; return hasher(ad.get_d());});

    pybind11::class_<ExactDouble> exact_double_class(module,"ExactDouble");
    exact_double_class.def(init<double>());
    exact_double_class.def(init<ExactDouble>());
    exact_double_class.def("__neg__", &__neg__<ExactDouble>);
    exact_double_class.def("__str__", &__cstr__<ExactDouble>);
    exact_double_class.def("__hash__", [](ExactDouble const& xd){std::hash<double> hasher; return hasher(xd.get_d());});

    module.def("exact", (ExactDouble(*)(double)) &cast_exact);
}


void export_integer(pymodule& module)
{
    pybind11::class_<Integer> integer_class(module,"Integer");
    integer_class.def(init<Int>());
    integer_class.def(init<Integer>());
    integer_class.def("__str__", &__cstr__<Integer>);
    integer_class.def("__repr__", &__repr__<Integer>);

    define_arithmetic(module,integer_class);
    define_lattice(module,integer_class);
    define_comparisons(module,integer_class);
    module.def("sqr", &_sqr_<Integer>);
    module.def("pow", &_pow_<Integer,Nat>);

    module.def("z", (Integer(*)(long long int)) &operator"" _z);

    implicitly_convertible<Int,Integer>();

    pybind11::class_<Natural,Integer> natural_class(module,"Natural");
    natural_class.def(init<Nat>());
    implicitly_convertible<Nat,Natural>();
}

void export_dyadic(pymodule& module)
{
    pybind11::class_<Dyadic> dyadic_class(module,"Dyadic");
    dyadic_class.def(init<Integer,Natural>());
    dyadic_class.def(init<ExactDouble>());
    dyadic_class.def(init<Integer>());
    dyadic_class.def(init<TwoExp>());
    dyadic_class.def(init<Dyadic>());
    dyadic_class.def(init<FloatDPValue>());
    dyadic_class.def(init<FloatMPValue>());
    dyadic_class.def("__str__", &__cstr__<Dyadic>);
    dyadic_class.def("__repr__", &__repr__<Dyadic>);
    dyadic_class.def("__rmul__", &__rmul__<Dyadic,Dyadic>);

    dyadic_class.def("__str__", &__cstr__<Dyadic>);
    define_infinitary(module,dyadic_class);
    define_arithmetic(module,dyadic_class);
    define_lattice(module,dyadic_class);
    define_comparisons(module,dyadic_class);
    module.def("sqr", &_sqr_<Dyadic>);
    module.def("hlf", &_hlf_<Dyadic>);

    module.def("dy", (Dyadic(*)(long double)) &operator""_dy);

    implicitly_convertible<Int,Dyadic>();
    implicitly_convertible<Integer,Dyadic>();
    implicitly_convertible<TwoExp,Dyadic>();
    implicitly_convertible<ExactDouble,Dyadic>();

    pybind11::class_<Two> two_class(module,"Two");
    two_class.def("__pow__", &__pow__<Two,Int>);
    two_class.def(__py_rdiv__, &__rdiv__<Two,Dyadic, Return<Dyadic>>);
    dyadic_class.def(__py_div__, &__div__<Dyadic,Two, Return<Dyadic>>);
    module.attr("two") = Two();

    pybind11::class_<TwoExp> two_exp_class(module,"TwoExp");
    two_exp_class.def(__py_rdiv__, &__rdiv__<TwoExp,Dyadic, Return<Dyadic>>);
    dyadic_class.def(__py_div__, &__div__<Dyadic,TwoExp, Return<Dyadic>>);
    two_exp_class.def("__str__", &__cstr__<Dyadic>);

    module.def("cast_exact", [](double d){return Dyadic(ExactDouble(d));});
}

void export_decimal(pymodule& module)
{
    pybind11::class_<Decimal> decimal_class(module,"Decimal");
    decimal_class.def(init<Decimal>());
    decimal_class.def(init<Dyadic>());
    decimal_class.def(init<Integer>());
    decimal_class.def(init<std::string>());
    decimal_class.def(init<double>());
    decimal_class.def("__str__", &__cstr__<Decimal>);
    decimal_class.def("__repr__", &__cstr__<Decimal>);
    define_arithmetic(module,decimal_class);
    define_lattice(module,decimal_class);
    define_comparisons(module,decimal_class);
    module.def("sqr", &_sqr_<Decimal>);
    module.def("hlf", &_hlf_<Decimal>);

    module.def("dec", (Decimal(*)(long double)) &operator"" _dec);

    implicitly_convertible<Int,Decimal>();
    implicitly_convertible<Integer,Decimal>();
    implicitly_convertible<Dyadic,Decimal>();
}

void export_rational(pymodule& module)
{
    pybind11::class_<Rational> rational_class(module,"Rational");
    rational_class.def(init<Integer,Integer>());
    rational_class.def(init<Int>());
    rational_class.def(init<Integer>());
    rational_class.def(init<Dyadic>());
    rational_class.def(init<Decimal>());
    rational_class.def(init<Rational>());
    rational_class.def("__str__", &__cstr__<Rational>);
    rational_class.def("__repr__", &__repr__<Rational>);

    define_infinitary(module,rational_class);

    define_arithmetic(module,rational_class);
    define_lattice(module,rational_class);
    define_comparisons(module,rational_class);
    module.def("hlf", &_hlf_<Rational>);
    module.def("sqr", &_sqr_<Rational>);
    module.def("rec", &_rec_<Rational>);

    module.def("q", (Rational(*)(long double)) &operator"" _q);

    implicitly_convertible<Int,Rational>();
    implicitly_convertible<Integer,Rational>();
    implicitly_convertible<Dyadic,Rational>();
    implicitly_convertible<Decimal,Rational>();
}

void export_real(pymodule& module)
{
    Real r; hlf(r);
    pybind11::class_<Real> real_class(module,"Real");
    real_class.def(init<Int>());
    real_class.def(init<Integer>());
    real_class.def(init<ExactDouble>());
    real_class.def(init<Dyadic>());
    real_class.def(init<Decimal>());
    real_class.def(init<Rational>());
    real_class.def(init<Real>());
    real_class.def("__str__", &__cstr__<Real>);
    real_class.def("__repr__", &__repr__<Real>);

    real_class.def("get", (FloatDPBounds(Real::*)(DoublePrecision)const) &Real::get);
    real_class.def("get", (FloatMPBounds(Real::*)(MultiplePrecision)const) &Real::get);
    real_class.def("compute", (ValidatedReal(Real::*)(Effort)const) &Real::compute);
    real_class.def("compute", (ValidatedReal(Real::*)(Accuracy)const) &Real::compute);

    define_elementary(module,real_class);
    define_lattice(module,real_class);
    define_comparisons(module,real_class);

    module.def("dist", &_dist_<Real,Real>);

    module.def("when", [](UpperKleenean p1, Real r1, UpperKleenean p2, Real r2){return when({p1,r1},{p2,r2});});
    module.def("choose", [](LowerKleenean p1, Real r1, LowerKleenean p2, Real r2){return choose({p1,r1},{p2,r2});});
    //module.def("when", (Real(*)(Case<UpperKleenean,Real>,Case<UpperKleenean,Real>)) &when);
    //module.def("choose", (Real(*)(Case<UpperKleenean,Real>,Case<UpperKleenean,Real>)) &choose);


    implicitly_convertible<Int,Real>();
    implicitly_convertible<Integer,Real>();
    implicitly_convertible<ExactDouble,Real>();
    implicitly_convertible<Dyadic,Real>();
    implicitly_convertible<Decimal,Real>();
    implicitly_convertible<Rational,Real>();

    pybind11::class_<PositiveReal,pybind11::bases<Real>> positive_real_class(module,"PositiveReal");
    positive_real_class.def("__repr__", &__cstr__<PositiveReal>);

    pybind11::class_<ValidatedReal> validated_real_class(module,"ValidatedReal");
    validated_real_class.def(init<DyadicBounds>());
    validated_real_class.def("__str__", &__cstr__<ValidatedReal>);
    validated_real_class.def("__repr__", &__cstr__<ValidatedReal>);

    validated_real_class.def("get", (DyadicBounds(ValidatedReal::*)()const) &ValidatedReal::get);
    validated_real_class.def("get", (FloatDPBounds(ValidatedReal::*)(DoublePrecision)const) &ValidatedReal::get);
    validated_real_class.def("get", (FloatMPBounds(ValidatedReal::*)(MultiplePrecision)const) &ValidatedReal::get);

    module.attr("pi") = pi;
}

template<class P> void export_number(pybind11::module& module)
{
    pybind11::class_<Number<P>> number_class(module,class_name<P>()+"Number");
    number_class.def(init<Rational>());
    number_class.def("class_name", &Number<P>::class_name);
    define_self_arithmetic(module,number_class);

//    number_class.def("get", (FloatDPBounds(Number<P>::*)(ValidatedTag,DoublePrecision)const) &Number<P>::get);
//    number_class.def("get", (FloatDPApproximation(Number<P>::*)(ApproximateTag,DoublePrecision)const) &Number<P>::get);
//    number_class.def("get", (FloatMPBounds(Number<P>::*)(ValidatedTag,MultiplePrecision)const) &Number<P>::get);
//    number_class.def("get", (FloatMPApproximation(Number<P>::*)(ApproximateTag,MultiplePrecision)const) &Number<P>::get);
}

FloatDPBounds get(ExactNumber const& y, DoublePrecision const& pr) { return y.get(pr); }
FloatMPBounds get(ExactNumber const& y, MultiplePrecision const& pr) { return y.get(pr); }
FloatDPBounds get(EffectiveNumber const& y, DoublePrecision const& pr) { return y.get(pr); }
FloatMPBounds get(EffectiveNumber const& y, MultiplePrecision const& pr) { return y.get(pr); }
FloatDPBounds get(ValidatedNumber const& y, DoublePrecision const& pr) { return y.get(pr); }
FloatMPBounds get(ValidatedNumber const& y, MultiplePrecision const& pr) { return y.get(pr); }
FloatDPUpperBound get(ValidatedUpperNumber const& y, DoublePrecision const& pr) { return y.get(pr); }
FloatMPUpperBound get(ValidatedUpperNumber const& y, MultiplePrecision const& pr) { return y.get(pr); }
FloatDPLowerBound get(ValidatedLowerNumber const& y, DoublePrecision const& pr) { return y.get(pr); }
FloatMPLowerBound get(ValidatedLowerNumber const& y, MultiplePrecision const& pr) { return y.get(pr); }
FloatDPApproximation get(ApproximateNumber const& y, DoublePrecision const& pr) { return y.get(pr); }
FloatMPApproximation get(ApproximateNumber const& y, MultiplePrecision const& pr) { return y.get(pr); }


template<class X, class Y> void implicitly_convertible_to() {
    if constexpr(Convertible<Y,ApproximateNumber>) { implicitly_convertible<X,ApproximateNumber>(); }
    if constexpr(Convertible<Y,ValidatedLowerNumber>) { implicitly_convertible<X,ValidatedLowerNumber>(); }
    if constexpr(Convertible<Y,ValidatedUpperNumber>) { implicitly_convertible<X,ValidatedUpperNumber>(); }
    if constexpr(Convertible<Y,ValidatedNumber>) { implicitly_convertible<X,ValidatedNumber>(); }
    if constexpr(Convertible<Y,EffectiveNumber>) { implicitly_convertible<X,EffectiveNumber>(); }
    if constexpr(Convertible<Y,ExactNumber>) { implicitly_convertible<X,ExactNumber>(); }
}

void export_numbers(pymodule& module)
{
    pybind11::class_<ExactNumber> exact_number_class(module,class_name<ExactNumber>().c_str());
    exact_number_class.def(init<ExactDouble>());
    exact_number_class.def(init<Rational>());
    exact_number_class.def(init<ExactNumber>());
    exact_number_class.def("get", (FloatDPBounds(*)(ExactNumber const&, DoublePrecision const&)) &get);
    exact_number_class.def("get", (FloatMPBounds(*)(ExactNumber const&, MultiplePrecision const&)) &get);
    exact_number_class.def("__str__", &__cstr__<ExactNumber>);
    exact_number_class.def("__repr__", &__repr__<ExactNumber>);


    pybind11::class_<EffectiveNumber> effective_number_class(module,class_name<EffectiveNumber>().c_str());
    effective_number_class.def(init<Rational>());
    effective_number_class.def(init<Real>());
    effective_number_class.def(init<ExactNumber>());
    effective_number_class.def(init<EffectiveNumber>());
    effective_number_class.def("get", (FloatDPBounds(*)(EffectiveNumber const&, DoublePrecision const&)) &get);
    effective_number_class.def("get", (FloatMPBounds(*)(EffectiveNumber const&, MultiplePrecision const&)) &get);
    effective_number_class.def("__str__", &__cstr__<EffectiveNumber>);
    effective_number_class.def("__repr__", &__repr__<EffectiveNumber>);

    define_elementary(module,effective_number_class);


    pybind11::class_<ValidatedNumber> validated_number_class(module,class_name<ValidatedNumber>().c_str());
    validated_number_class.def(init<Rational>());
    validated_number_class.def(init<Real>());
    validated_number_class.def(init<DyadicBounds>());
    validated_number_class.def(init<ExactNumber>());
    validated_number_class.def(init<EffectiveNumber>());
    validated_number_class.def(init<ValidatedNumber>());
    validated_number_class.def("get", (FloatDPBounds(*)(ValidatedNumber const&, DoublePrecision const&)) &get);
    validated_number_class.def("get", (FloatMPBounds(*)(ValidatedNumber const&, MultiplePrecision const&)) &get);
    validated_number_class.def("__str__", &__cstr__<ValidatedNumber>);
    validated_number_class.def("__repr__", &__repr__<ValidatedNumber>);

    define_elementary(module,validated_number_class);
    validated_number_class.def(pybind11::init([](pybind11::dict pydct){return ValidatedNumber(bounds_from_dict<DyadicBounds>(pydct));}));
    implicitly_convertible<pybind11::dict,ValidatedNumber>();


    pybind11::class_<ValidatedUpperNumber> validated_upper_number_class(module,class_name<ValidatedUpperNumber>().c_str());
    validated_upper_number_class.def(init<ValidatedNumber>());
    validated_upper_number_class.def("get", (FloatDPUpperBound(*)(ValidatedUpperNumber const&, DoublePrecision const&)) &get);
    validated_upper_number_class.def("get", (FloatMPUpperBound(*)(ValidatedUpperNumber const&, MultiplePrecision const&)) &get);
    validated_upper_number_class.def("__str__", &__cstr__<ValidatedUpperNumber>);
    validated_upper_number_class.def("__repr__", &__repr__<ValidatedUpperNumber>);

    define_monotonic(module,validated_upper_number_class);


    pybind11::class_<ValidatedLowerNumber> validated_lower_number_class(module,(class_name<ValidatedLowerNumber>()).c_str());
    validated_lower_number_class.def(init<ValidatedNumber>());
    validated_lower_number_class.def("get", (FloatDPLowerBound(*)(ValidatedLowerNumber const&, DoublePrecision const&)) &get);
    validated_lower_number_class.def("get", (FloatMPLowerBound(*)(ValidatedLowerNumber const&, MultiplePrecision const&)) &get);
    validated_lower_number_class.def("__str__", &__cstr__<ValidatedLowerNumber>);
    validated_lower_number_class.def("__repr__", &__repr__<ValidatedLowerNumber>);

    define_monotonic(module,validated_lower_number_class);


    pybind11::class_<ApproximateNumber> approximate_number_class(module,class_name<ApproximateNumber>().c_str());
    approximate_number_class.def(init<Int>());
    approximate_number_class.def(init<Dbl>());
    approximate_number_class.def(init<Rational>());
    approximate_number_class.def(init<Real>());
    approximate_number_class.def(init<ExactNumber>());
    approximate_number_class.def(init<EffectiveNumber>());
    approximate_number_class.def(init<ValidatedNumber>());
    approximate_number_class.def(init<ApproximateNumber>());
    approximate_number_class.def("get", (FloatDPApproximation(*)(ApproximateNumber const&, DoublePrecision const&)) &get);
    approximate_number_class.def("get", (FloatMPApproximation(*)(ApproximateNumber const&, MultiplePrecision const&)) &get);
    approximate_number_class.def("__str__", &__cstr__<ApproximateNumber>);
    approximate_number_class.def("__repr__", &__repr__<ApproximateNumber>);

    define_elementary(module,approximate_number_class);
    define_lattice(module,approximate_number_class);


    implicitly_convertible<ValidatedNumber,ApproximateNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedLowerNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedUpperNumber>();
    implicitly_convertible<EffectiveNumber,ApproximateNumber>();
    implicitly_convertible<EffectiveNumber,ValidatedLowerNumber>();
    implicitly_convertible<EffectiveNumber,ValidatedUpperNumber>();
    implicitly_convertible<EffectiveNumber,ValidatedNumber>();
    implicitly_convertible<ExactNumber,ApproximateNumber>();
    implicitly_convertible<ExactNumber,ValidatedLowerNumber>();
    implicitly_convertible<ExactNumber,ValidatedUpperNumber>();
    implicitly_convertible<ExactNumber,ValidatedNumber>();
    implicitly_convertible<ExactNumber,EffectiveNumber>();

    implicitly_convertible<DyadicBounds,ValidatedNumber>();

    implicitly_convertible_to<Real,EffectiveNumber>();
    implicitly_convertible_to<Rational,ExactNumber>();
    implicitly_convertible_to<Decimal,ExactNumber>();
    implicitly_convertible_to<Dyadic,ExactNumber>();
    implicitly_convertible_to<ExactDouble,ExactNumber>();
    implicitly_convertible_to<Integer,ExactNumber>();
    implicitly_convertible_to<Int,ExactNumber>();

    implicitly_convertible_to<double,ApproximateNumber>();


    // TODO: These exports should be with FloatXxx
    exact_number_class.def(init<FloatDPValue>());
    exact_number_class.def(init<FloatMPValue>());
    validated_number_class.def(init<FloatDPBall>());
    validated_number_class.def(init<FloatMPDPBall>());
    validated_number_class.def(init<FloatMPBall>());
    validated_number_class.def(init<FloatDPBounds>());
    validated_number_class.def(init<FloatMPBounds>());
    validated_upper_number_class.def(init<FloatDPUpperBound>());
    validated_upper_number_class.def(init<FloatMPUpperBound>());
    validated_lower_number_class.def(init<FloatDPLowerBound>());
    validated_lower_number_class.def(init<FloatMPLowerBound>());
    approximate_number_class.def(init<FloatDPApproximation>());
    approximate_number_class.def(init<FloatMPApproximation>());

    define_mixed_arithmetic<ValidatedNumber,FloatBall<DP,DP>>(module,validated_number_class);
    define_mixed_arithmetic<ValidatedNumber,FloatBall<MP,DP>>(module,validated_number_class);
    define_mixed_arithmetic<ValidatedNumber,FloatBall<MP,MP>>(module,validated_number_class);
    define_mixed_arithmetic<ValidatedNumber,FloatBounds<DP>>(module,validated_number_class);
    define_mixed_arithmetic<ValidatedNumber,FloatBounds<MP>>(module,validated_number_class);
    define_mixed_arithmetic<ApproximateNumber,FloatApproximation<DP>>(module,approximate_number_class);
    define_mixed_arithmetic<ApproximateNumber,FloatApproximation<MP>>(module,approximate_number_class);
}




void export_dyadic_bounds(pymodule& module)
{
    pybind11::class_<DyadicBounds> dyadic_bounds_class(module,"DyadicBounds");
    dyadic_bounds_class.def(init<DyadicBounds>());
    dyadic_bounds_class.def(init<Dyadic>());
    dyadic_bounds_class.def(init<Dyadic,Dyadic>());
    dyadic_bounds_class.def(pybind11::init([](pybind11::dict pydct){return bounds_from_dict<DyadicBounds>(pydct);}));
    implicitly_convertible<int,DyadicBounds>();
    implicitly_convertible<Dyadic,DyadicBounds>();
    implicitly_convertible<pybind11::dict,DyadicBounds>();
    dyadic_bounds_class.def("lower", &DyadicBounds::lower);
    dyadic_bounds_class.def("upper", &DyadicBounds::upper);

    define_arithmetic(module,dyadic_bounds_class);
    define_lattice(module,dyadic_bounds_class);
    dyadic_bounds_class.def("hlf", &_hlf_<DyadicBounds>);

    dyadic_bounds_class.def("__str__", &__cstr__<DyadicBounds>);
    dyadic_bounds_class.def("__repr__", &__repr__<DyadicBounds>);
}

void export_rational_bounds(pymodule& module)
{
    pybind11::class_<RationalBounds> rational_bounds_class(module,"RationalBounds");
    rational_bounds_class.def(init<RationalBounds>());
    rational_bounds_class.def(init<Rational>());
    rational_bounds_class.def(init<Rational,Rational>());
    rational_bounds_class.def(init<DyadicBounds>());
    rational_bounds_class.def(pybind11::init([](pybind11::dict pydct){return bounds_from_dict<RationalBounds>(pydct);}));
    implicitly_convertible<int,RationalBounds>();
    implicitly_convertible<Rational,RationalBounds>();
    implicitly_convertible<pybind11::dict,RationalBounds>();
    rational_bounds_class.def("lower", &RationalBounds::lower);
    rational_bounds_class.def("upper", &RationalBounds::upper);

    define_arithmetic(module,rational_bounds_class);
    define_lattice(module,rational_bounds_class);
    module.def("sqr", &_sqr_<RationalBounds>);
    module.def("rec", &_rec_<RationalBounds>);

    rational_bounds_class.def("__str__", &__cstr__<RationalBounds>);
    rational_bounds_class.def("__repr__", &__repr__<RationalBounds>);
}


namespace Ariadne {
const Rounding round_upward = Rounding(upward);
const Rounding round_downward = Rounding(downward);
const Rounding round_to_nearest = Rounding(to_nearest);
}

void export_rounding_mode(pymodule& module) {
    pybind11::class_<Rounding> rounding_mode_class(module,"Rounding");
    rounding_mode_class.def(init<Rounding>());
    rounding_mode_class.def("__str__", &__cstr__<Rounding>);
    module.attr("upward") = round_upward;
    module.attr("downward") = round_downward;
    module.attr("to_near") = round_to_nearest;
    module.attr("up") = round_upward;
    module.attr("down") = round_downward;
    module.attr("near") = round_to_nearest;
}


template<class PR> Void export_precision(pymodule& module);

template<> Void export_precision<DoublePrecision>(pymodule& module) {
    pybind11::class_<DoublePrecision> precision_class(module,"DoublePrecision");
    precision_class.def(init<>());
    precision_class.def("__str__", &__cstr__<DoublePrecision>);
    precision_class.def("__repr__", &__cstr__<DoublePrecision>);
    module.attr("double_precision") = double_precision;
    module.attr("dp") = dp;
}

template<> Void export_precision<MultiplePrecision>(pymodule& module) {
    pybind11::class_<MultiplePrecision> precision_class(module,"MultiplePrecision");
    precision_class.def(init<Nat>());
    precision_class.def("bits",&MultiplePrecision::bits);
    precision_class.def("__str__", &__cstr__<MultiplePrecision>);
    precision_class.def("__repr__", &__cstr__<MultiplePrecision>);
    module.def("multiple_precision", (MultiplePrecision(*)(mpfr_prec_t)) &precision);
    module.def("precision", (MultiplePrecision(*)(mpfr_prec_t)) &precision);
}


template<class PR> void export_raw_float(pymodule& module)
{
    typedef RawFloat<PR> F;

//    typedef typename F::RoundingModeType RND;
//    implicitly_convertible<Rounding,RND>();
    typedef Rounding RND;

    pybind11::class_<F> raw_float_class(module,("Float"+numeric_class_tag<PR>()).c_str());
    raw_float_class.def(pybind11::init([](double d, PR pr){return F(cast_exact(d),pr); }));
    raw_float_class.def(init<ExactDouble,PR>());
    raw_float_class.def(init<Dyadic,PR>());
    raw_float_class.def(init<Rational,RND,PR>());

    define_infinitary(module,raw_float_class);

    raw_float_class.def_static("eps", (F(*)(PR)) &F::eps);
    raw_float_class.def_static("max", (F(*)(PR)) &F::max);
    raw_float_class.def_static("min", (F(*)(PR)) &F::min);

    raw_float_class.def("__str__", &__cstr__<RawFloat<PR>>);
    raw_float_class.def("__repr__", &__repr__<RawFloat<PR>>);

    module.def("nul", &_nul_<F>);
    module.def("pos", &_pos_<F>);
    module.def("neg", &_neg_<F>);
    module.def("hlf", &_hlf_<F>);

    module.def("nul", &_nul_<RND,F>);
    module.def("pos", &_pos_<RND,F>);
    module.def("neg", &_neg_<RND,F>);
    module.def("hlf", &_hlf_<RND,F>);
    module.def("sqr", &_sqr_<RND,F>);
    module.def("rec", &_rec_<RND,F>);
    module.def("add", &_add_<RND,F,F>);
    module.def("sub", &_sub_<RND,F,F>);
    module.def("mul", &_mul_<RND,F,F>);
    module.def("div", &_div_<RND,F,F>);
    module.def("fma", &_fma_<RND,F,F,F>);
    module.def("pow", &_pow_<RND,F,Int>);
    module.def("sqrt", &_sqrt_<RND,F>);
    module.def("exp", &_exp_<RND,F>);
    module.def("log", &_log_<RND,F>);
    module.def("sin", &_sin_<RND,F>);
    module.def("cos", &_cos_<RND,F>);
    module.def("tan", &_tan_<RND,F>);
    module.def("asin", &_asin_<RND,F>);
    module.def("acos", &_acos_<RND,F>);
    module.def("atan", &_atan_<RND,F>);

    module.def("abs", &_abs_<F>);
    module.def("max", &_min_<F,F>);
    module.def("min", &_max_<F,F>);

    define_comparisons<F>(module,raw_float_class);

}



template<class PR> void export_float_value(pymodule& module)
{
    pybind11::class_<FloatValue<PR>> float_value_class(module,("Float"+numeric_class_tag<PR>()+"Value").c_str());
    float_value_class.def(init<PR>());
    float_value_class.def(init<RawFloat<PR>>());
    float_value_class.def(init<ExactDouble,PR>());
    float_value_class.def(init<Integer,PR>());
    float_value_class.def(init<Dyadic,PR>());
    float_value_class.def(init<FloatValue<PR>>());

    float_value_class.def("__pos__", &__pos__<FloatValue<PR>>);
    float_value_class.def("__neg__", &__neg__<FloatValue<PR>>);

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    define_mixed_arithmetic<FloatValue<PR>,Int>(module,float_value_class);
    define_mixed_arithmetic<FloatValue<PR>,Dyadic>(module,float_value_class);
    define_elementary<FloatValue<PR>>(module,float_value_class);
    //float_value_class.define_mixed_arithmetic<ApproximateNumericType>();
    //float_value_class.define_mixed_arithmetic<LowerNumericType>();
    //float_value_class.define_mixed_arithmetic<UpperNumericType>();
    //float_value_class.define_mixed_arithmetic<ValidatedNumericType>();
    define_lattice(module,float_value_class);
    define_mixed_lattice<FloatValue<PR>,Int>(module,float_value_class);
    define_mixed_lattice<FloatValue<PR>,Dyadic>(module,float_value_class);
    define_comparisons(module,float_value_class);

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("raw", (RawFloat<PR>const&(FloatValue<PR>::*)()const)&FloatValue<PR>::raw);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    float_value_class.def("__str__", &__cstr__<FloatValue<PR>>);
    float_value_class.def("__repr__", &__repr__<FloatValue<PR>>);

    float_value_class.def_static("set_output_places",&FloatValue<PR>::set_output_places);
}


template<class PRE> void export_float_error(pymodule& module)
{
    pybind11::class_<FloatError<PRE>> float_error_class(module,("Float"+numeric_class_tag<PRE>()+"Error").c_str());
    float_error_class.def(init<PRE>());
    float_error_class.def(init<RawFloat<PRE>>());
    float_error_class.def(init<Nat,PRE>());
    float_error_class.def(init<FloatError<PRE>>());

    float_error_class.def("__pos__", &__pos__<FloatError<PRE>>, pybind11::is_operator());
    float_error_class.def("__add__", &__add__<FloatError<PRE>,FloatError<PRE>>, pybind11::is_operator());
    float_error_class.def("__mul__", &__mul__<FloatError<PRE>,FloatError<PRE>>, pybind11::is_operator());

    float_error_class.def("__add__", &__add__<FloatError<PRE>,Nat>, pybind11::is_operator());
    float_error_class.def("__radd__", &__radd__<FloatError<PRE>,Nat>, pybind11::is_operator());
    float_error_class.def("__rmul__", &__rmul__<FloatError<PRE>,Nat>, pybind11::is_operator());
    float_error_class.def("__mul__", &__mul__<FloatError<PRE>,Nat>, pybind11::is_operator());
    float_error_class.def(__py_div__, &__div__<FloatError<PRE>,Nat>, pybind11::is_operator());

    module.def("max", &_max_<FloatError<PRE>,FloatError<PRE>>);
    module.def("min", &_min_<FloatError<PRE>,FloatError<PRE>>);

    float_error_class.def("raw",(RawFloat<PRE>const&(FloatError<PRE>::*)()const)&FloatError<PRE>::raw);
    float_error_class.def("__str__", &__cstr__<FloatError<PRE>>);
    float_error_class.def("__repr__", &__repr__<FloatError<PRE>>);

    float_error_class.def("precision", &FloatError<PRE>::precision);

    module.def("log2", (FloatUpperBound<PRE>(*)(FloatError<PRE>const&)) &_log2_);
    float_error_class.def_static("set_output_places",&FloatError<PRE>::set_output_places);
}


template<class PR, class PRE=PR> void export_float_ball(pymodule& module)
{
    String tag = Same<PR,PRE> ? numeric_class_tag<PR>() : numeric_class_tag<PR>()+numeric_class_tag<PRE>();
    pybind11::class_<FloatBall<PR,PRE>> float_ball_class(module,("Float"+tag+"Ball").c_str());
    float_ball_class.def(init<PR,PRE>());
    float_ball_class.def(init<RawFloat<PR>,RawFloat<PRE>>());
    float_ball_class.def(init<FloatValue<PR>,FloatError<PRE>>());
    float_ball_class.def(init<Real,PR>());

    float_ball_class.def(init<ExactDouble,PR>());
    float_ball_class.def(init<ValidatedNumber,PR>());
    float_ball_class.def(init<FloatValue<PR>>());
    float_ball_class.def(init<FloatBall<PR,PRE>>());
    float_ball_class.def(init<FloatBounds<PR>>());

    float_ball_class.def("value", &FloatBall<PR,PRE>::value);
    float_ball_class.def("error", &FloatBall<PR,PRE>::error);
    float_ball_class.def("lower", &FloatBall<PR,PRE>::lower);
    float_ball_class.def("upper", &FloatBall<PR,PRE>::upper);

    float_ball_class.def("precision", &FloatBall<PR,PRE>::precision);
    float_ball_class.def("error_precision", &FloatBall<PR,PRE>::error_precision);

    static_assert(Same<decltype(declval<FloatBall<PR,PRE>>()+declval<ValidatedNumber>()),FloatBall<PR,PRE>>);
    static_assert(Same<decltype(add(declval<FloatBall<PR,PRE>>(),declval<ValidatedNumber>())),FloatBall<PR,PRE>>);
    static_assert(Same<decltype(max(declval<FloatBall<PR,PRE>>(),declval<ValidatedNumber>())),FloatBall<PR,PRE>>);

    define_mixed_arithmetic<FloatBall<PR,PRE>,ValidatedNumber>(module,float_ball_class);
    define_mixed_arithmetic<FloatBall<PR,PRE>,Int>(module,float_ball_class);
    define_elementary<FloatBall<PR,PRE>>(module,float_ball_class);
    define_lattice<FloatBall<PR,PRE>>(module,float_ball_class);
    define_mixed_lattice<FloatBall<PR,PRE>,ValidatedNumber>(module,float_ball_class);
    define_mixed_lattice<FloatBall<PR,PRE>,Int>(module,float_ball_class);
    define_comparisons<FloatBall<PR,PRE>>(module,float_ball_class);
    define_mixed_comparisons<FloatBall<PR,PRE>,ValidatedNumber>(module,float_ball_class);
//    float_ball_class.define_mixed_arithmetic<FloatBall<PR,PRE>>();
//    float_ball_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_ball_class.define_mixed_arithmetic<LowerNumericType>();
//    float_ball_class.define_mixed_arithmetic<UpperNumericType>();
//    float_ball_class.define_mixed_arithmetic<ValidatedNumericType>();

    define_mixed_arithmetic<FloatBall<PR,PRE>,FloatBounds<PR>>(module,float_ball_class);
    define_mixed_lattice<FloatBall<PR,PRE>,FloatBounds<PR>>(module,float_ball_class);

    float_ball_class.def("__str__", &__cstr__<FloatBall<PR,PRE>>);
    float_ball_class.def("__repr__", &__repr__<FloatBall<PR,PRE>>);
}


template<class PR> void export_float_bounds(pymodule& module)
{
    pybind11::class_<FloatBounds<PR>> float_bounds_class(module,("Float"+numeric_class_tag<PR>()+"Bounds").c_str());
    float_bounds_class.def(init<PR>());
    float_bounds_class.def(init<RawFloat<PR>>());
    float_bounds_class.def(init<RawFloat<PR>,RawFloat<PR>>());
    float_bounds_class.def(init<FloatLowerBound<PR>,FloatUpperBound<PR>>());
    float_bounds_class.def(init<ValidatedLowerNumber,ValidatedUpperNumber,PR>());
    float_bounds_class.def(init<Real,PR>());

    float_bounds_class.def(init<ExactDouble,PR>());
    float_bounds_class.def(init<ExactDouble,ExactDouble,PR>());
    float_bounds_class.def(init<ValidatedNumber,PR>());
    float_bounds_class.def(init<FloatValue<PR>>());
    float_bounds_class.def(init<FloatBall<PR>>());
    float_bounds_class.def(init<FloatBounds<PR>>());

    float_bounds_class.def("__str__", &__cstr__<FloatBounds<PR>>);
    float_bounds_class.def("__repr__", &__repr__<FloatBounds<PR>>);

    float_bounds_class.def("precision", &FloatBounds<PR>::precision);
    float_bounds_class.def("lower", &FloatBounds<PR>::lower);
    float_bounds_class.def("upper", &FloatBounds<PR>::upper);
    float_bounds_class.def("value", &FloatBounds<PR>::value);
    float_bounds_class.def("error", &FloatBounds<PR>::error);

    static_assert(Same<decltype(declval<FloatBounds<PR>>()+declval<ValidatedNumber>()),FloatBounds<PR>>);
    static_assert(Same<decltype(add(declval<FloatBounds<PR>>(),declval<ValidatedNumber>())),FloatBounds<PR>>);
    static_assert(Same<decltype(max(declval<FloatBounds<PR>>(),declval<ValidatedNumber>())),FloatBounds<PR>>);
    static_assert(Same<decltype(declval<ValidatedNumber>()+declval<FloatBounds<PR>>()),FloatBounds<PR>>);
    static_assert(Same<decltype(add(declval<ValidatedNumber>(),declval<FloatBounds<PR>>())),FloatBounds<PR>>);
    static_assert(Same<decltype(max(declval<ValidatedNumber>(),declval<FloatBounds<PR>>())),FloatBounds<PR>>);

    static_assert(Same<decltype(declval<FloatBounds<PR>>()+declval<FloatBall<PR>>()),FloatBounds<PR>>);
    static_assert(Same<decltype(add(declval<FloatBounds<PR>>(),declval<FloatBall<PR>>())),FloatBounds<PR>>);
    static_assert(Same<decltype(max(declval<FloatBounds<PR>>(),declval<FloatBall<PR>>())),FloatBounds<PR>>);
    static_assert(Same<decltype(declval<FloatBall<PR>>()+declval<FloatBounds<PR>>()),FloatBounds<PR>>);
    static_assert(Same<decltype(add(declval<FloatBall<PR>>(),declval<FloatBounds<PR>>())),FloatBounds<PR>>);
    static_assert(Same<decltype(max(declval<FloatBall<PR>>(),declval<FloatBounds<PR>>())),FloatBounds<PR>>);

    define_mixed_arithmetic<FloatBounds<PR>,ValidatedNumber>(module,float_bounds_class);
    define_mixed_arithmetic<FloatBounds<PR>,Int>(module,float_bounds_class);
    define_mixed_arithmetic<FloatBounds<PR>,FloatBall<PR>>(module,float_bounds_class);
    define_elementary<FloatBounds<PR>>(module,float_bounds_class);
    define_lattice<FloatBounds<PR>>(module,float_bounds_class);
    define_mixed_lattice<FloatBounds<PR>,ValidatedNumber>(module,float_bounds_class);
    define_mixed_lattice<FloatBounds<PR>,Int>(module,float_bounds_class);
    define_comparisons<FloatBounds<PR>>(module,float_bounds_class);
    define_mixed_comparisons<FloatBounds<PR>,ValidatedNumber>(module,float_bounds_class);
//    float_bounds_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_bounds_class.define_mixed_arithmetic<LowerNumericType>();
//    float_bounds_class.define_mixed_arithmetic<UpperNumericType>();
//    float_bounds_class.define_mixed_arithmetic<ValidatedNumericType>();

    module.def("refinement", &_refinement_<FloatBounds<PR>,FloatBounds<PR>>);
    module.def("refines", &_refines_<FloatBounds<PR>,FloatBounds<PR>>);
    module.def("inconsistent", &_inconsistent_<FloatBounds<PR>,FloatBounds<PR>>);

    implicitly_convertible<FloatValue<PR>,FloatBounds<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatBounds<PR>>();

    float_bounds_class.def_static("set_output_places",&FloatBounds<PR>::set_output_places);
}


template<class PR> void export_float_upper_bound(pymodule& module)
{
    pybind11::class_<FloatUpperBound<PR>> float_upper_bound_class(module,("Float"+numeric_class_tag<PR>()+"UpperBound").c_str());
    float_upper_bound_class.def(init<PR>());
    float_upper_bound_class.def(init<RawFloat<PR>>());
    float_upper_bound_class.def(init<ExactDouble,PR>());
    float_upper_bound_class.def(init<ValidatedUpperNumber,PR>());

    float_upper_bound_class.def(init<ExactDouble,PR>());
    float_upper_bound_class.def(init<FloatValue<PR>>());
    float_upper_bound_class.def(init<FloatBall<PR>>());
    float_upper_bound_class.def(init<FloatBounds<PR>>());
    float_upper_bound_class.def(init<FloatUpperBound<PR>>());
    float_upper_bound_class.def(init<FloatError<PR>>());
    float_upper_bound_class.def("__str__", &__cstr__<FloatUpperBound<PR>>);
    float_upper_bound_class.def("__repr__", &__repr__<FloatUpperBound<PR>>);

    float_upper_bound_class.def("precision", &FloatUpperBound<PR>::precision);
    float_upper_bound_class.def("raw", (RawFloat<PR>const&(FloatUpperBound<PR>::*)()const)&FloatUpperBound<PR>::raw);

    define_monotonic(module,float_upper_bound_class);
    define_mixed_monotonic<FloatUpperBound<PR>,Int>(module,float_upper_bound_class);
    define_mixed_monotonic<FloatUpperBound<PR>,ValidatedUpperNumber>(module,float_upper_bound_class);
    define_lattice<FloatUpperBound<PR>>(module,float_upper_bound_class);
    define_mixed_lattice<FloatUpperBound<PR>,ValidatedUpperNumber>(module,float_upper_bound_class);

    float_upper_bound_class.def("__lt__", &__lt__<FloatUpperBound<PR>,FloatLowerBound<PR>>);
    float_upper_bound_class.def("__le__", &__le__<FloatUpperBound<PR>,FloatLowerBound<PR>>);

//    float_upper_bound_class.def(self > FloatBounds<PR>());
//    float_upper_bound_class.def(self > FloatApproximation<PR>());
//    float_upper_bound_class.def(self < FloatLowerBound<PR>());
//    float_upper_bound_class.def(self >= FloatBounds<PR>());
//    float_upper_bound_class.def(self >= FloatApproximation<PR>());
//    float_upper_bound_class.def(self <= FloatLowerBound<PR>());

//    float_upper_bound_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_upper_bound_class.def(UpperNumericType() + self);
//    float_upper_bound_class.def(LowerNumericType() - self);
//    float_upper_bound_class.def(self + UpperNumericType());
//    float_upper_bound_class.def(self - LowerNumericType());

    implicitly_convertible<FloatBounds<PR>,FloatUpperBound<PR>>();
    implicitly_convertible<FloatError<PR>,FloatUpperBound<PR>>();
}


template<class PR> void export_float_lower_bound(pymodule& module)
{
    pybind11::class_<FloatLowerBound<PR>> float_lower_bound_class(module,("Float"+numeric_class_tag<PR>()+"LowerBound").c_str());
    float_lower_bound_class.def(init<PR>());
    float_lower_bound_class.def(init<RawFloat<PR>>());
    float_lower_bound_class.def(init<ExactDouble,PR>());
    float_lower_bound_class.def(init<ValidatedLowerNumber,PR>());

    float_lower_bound_class.def(init<FloatValue<PR>>());
    float_lower_bound_class.def(init<FloatBall<PR>>());
    float_lower_bound_class.def(init<FloatBounds<PR>>());
    float_lower_bound_class.def(init<FloatLowerBound<PR>>());
    float_lower_bound_class.def("__str__", &__cstr__<FloatLowerBound<PR>>);
    float_lower_bound_class.def("__repr__", &__repr__<FloatLowerBound<PR>>);

    float_lower_bound_class.def("precision", &FloatLowerBound<PR>::precision);

    define_monotonic(module,float_lower_bound_class);
    define_mixed_monotonic<FloatLowerBound<PR>,Int>(module,float_lower_bound_class);
    define_mixed_monotonic<FloatLowerBound<PR>,ValidatedLowerNumber>(module,float_lower_bound_class);
    define_lattice<FloatLowerBound<PR>>(module,float_lower_bound_class);
    define_mixed_lattice<FloatLowerBound<PR>,ValidatedLowerNumber>(module,float_lower_bound_class);

    float_lower_bound_class.def("__lt__", &__lt__<FloatLowerBound<PR>,FloatUpperBound<PR>>);
    float_lower_bound_class.def("__le__", &__le__<FloatLowerBound<PR>,FloatUpperBound<PR>>);

//    float_lower_bound_class.def(self < FloatBounds<PR>());
//    float_lower_bound_class.def(self < FloatApproximation<PR>());
//    float_lower_bound_class.def(self > FloatUpperBound<PR>());
//    float_lower_bound_class.def(self <= FloatBounds<PR>());
//    float_lower_bound_class.def(self <= FloatApproximation<PR>());
//    float_lower_bound_class.def(self >= FloatUpperBound<PR>());
//    float_lower_bound_class.def("raw", (RawFloat<PR>const&(FloatLowerBound<PR>::*)()const)&FloatLowerBound<PR>::raw);

//    float_lower_bound_class.define_mixed_arithmetic<ApproximateNumericType>();

//    float_lower_bound_class.def(LowerNumericType() + self);
//    float_lower_bound_class.def(UpperNumericType() - self);
//    float_lower_bound_class.def(self + LowerNumericType());
//    float_lower_bound_class.def(self - UpperNumericType());

    implicitly_convertible<FloatBounds<PR>,FloatLowerBound<PR>>();
}


template<class PR> void export_float_approximation(pymodule& module)
{
    pybind11::class_<FloatApproximation<PR>> float_approximation_class(module,("Float"+numeric_class_tag<PR>()+"Approximation").c_str());

    float_approximation_class.def(init<PR>());
    float_approximation_class.def(init<RawFloat<PR>>());
    float_approximation_class.def(init<double,PR>());
    float_approximation_class.def(init<ExactDouble,PR>());
    float_approximation_class.def(init<Real,PR>());
    float_approximation_class.def(init<ApproximateNumber,PR>());

    float_approximation_class.def(init<FloatValue<PR>>());
    float_approximation_class.def(init<FloatBall<PR>>());
    float_approximation_class.def(init<FloatBounds<PR>>());
    float_approximation_class.def(init<FloatLowerBound<PR>>());
    float_approximation_class.def(init<FloatUpperBound<PR>>());
    float_approximation_class.def(init<FloatApproximation<PR>>());

    float_approximation_class.def("precision", &FloatApproximation<PR>::precision);
    float_approximation_class.def("raw", (RawFloat<PR>const&(FloatApproximation<PR>::*)()const)&FloatApproximation<PR>::raw);
    float_approximation_class.def("get_d", &FloatApproximation<PR>::get_d);

    define_mixed_arithmetic<FloatApproximation<PR>,double>(module,float_approximation_class);
    define_mixed_arithmetic<FloatApproximation<PR>,ApproximateNumber>(module,float_approximation_class);
    define_elementary<FloatApproximation<PR>>(module,float_approximation_class);
    define_lattice<FloatApproximation<PR>>(module,float_approximation_class);
    define_mixed_lattice<FloatApproximation<PR>,ApproximateNumber>(module,float_approximation_class);
    define_comparisons<FloatApproximation<PR>>(module,float_approximation_class);
    define_mixed_comparisons<FloatApproximation<PR>,ApproximateNumber>(module,float_approximation_class);

    float_approximation_class.def("__str__", &__cstr__<FloatApproximation<PR>>);
    float_approximation_class.def("__repr__", &__repr__<FloatApproximation<PR>>);

    float_approximation_class.def_static("set_output_places",&FloatApproximation<PR>::set_output_places);

    //NOTE: Conversion to FloatMP needs precision, so disallow conversion to FloatDP as well
    //implicitly_convertible<double,FloatApproximation<PR>>();
    implicitly_convertible<FloatValue<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBounds<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatUpperBound<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatLowerBound<PR>,FloatApproximation<PR>>();

    module.def("cast_exact",(FloatValue<PR>const&(*)(FloatApproximation<PR> const&)) &cast_exact);
    module.def("cast_exact",(FloatValue<PR>const&(*)(RawFloat<PR> const&)) &cast_exact);
}


template<class PR> Void export_user_floats(pymodule& module) {
    export_float_error<PR>(module);
    export_float_value<PR>(module);
    export_float_ball<PR>(module);
    export_float_bounds<PR>(module);
    export_float_upper_bound<PR>(module);
    export_float_lower_bound<PR>(module);
    export_float_approximation<PR>(module);

    module.def("Error", [](Nat m, PR pr){return Error(m,pr);});
    module.def("Error", [](RawFloatType<PR> const& v){return Error(v);});

    module.def("Value", [](Dyadic const& y, PR pr){return Value(y,pr);});
    module.def("Value", [](RawFloatType<PR> const& v){return Value(v);});

    if constexpr (Same<PR,MultiplePrecision>) {
        module.def("Ball", [](EffectiveNumber const& y, MultiplePrecision pr, DoublePrecision pre){return Ball(y,pr,pre);});
        module.def("Ball", [](ValidatedNumber const& y, MultiplePrecision pr, DoublePrecision pre){return Ball(y,pr,pre);});
        module.def("Ball", [](RawFloatType<MultiplePrecision> const& v, RawFloatType<DoublePrecision> const& e){return Ball(v,e);});
    }

    module.def("Ball", [](EffectiveNumber const& y, PR pr){return Ball<RawFloatType<PR>>(y,pr);});
    module.def("Ball", [](ValidatedNumber const& y, PR pr){return Ball(y,pr);});
    module.def("Ball", [](ValidatedNumber const& y, PR pr, PR pre){return Ball(y,pr,pre);});
    module.def("Ball", [](RawFloatType<PR> const& v, RawFloatType<PR> const& e){return Ball(v,e);});
    module.def("Bounds", [](EffectiveNumber const& y, PR pr){return Bounds(y,pr);});
    module.def("Bounds", [](ValidatedNumber const& y, PR pr){return Bounds(y,pr);});
    module.def("Bounds", [](ValidatedLowerNumber const& yl, ValidatedUpperNumber const& yu, PR pr){return Bounds(yl,yu,pr);});
    module.def("Bounds", [](RawFloatType<PR> const& l, RawFloatType<PR> const& u){return Bounds(l,u);});
    module.def("UpperBound", [](ValidatedUpperNumber const& y, PR pr){return UpperBound(y,pr);});
    module.def("UpperBound", [](RawFloatType<PR> const& u){return UpperBound(u);});
    module.def("LowerBound", [](ValidatedLowerNumber const& y, PR pr){return LowerBound(y,pr);});
    module.def("LowerBound", [](RawFloatType<PR> const& l){return LowerBound(l);});
    module.def("Approximation", [](ApproximateNumber const& y, PR pr){return Approximation(y,pr);});
    module.def("Approximation", [](RawFloatType<PR> const& a){return Approximation(a);});
}


template<class PR> Void export_float_to_number_conversions(pymodule& module) {
    if constexpr(ALLOW_CONCRETE_TO_GENERIC_NUMBER_CONVERSIONS) {
        implicitly_convertible<FloatValue<PR>,ExactNumber>();
        implicitly_convertible<FloatValue<PR>,ValidatedNumber>();
        implicitly_convertible<FloatValue<PR>,ValidatedUpperNumber>();
        implicitly_convertible<FloatValue<PR>,ValidatedLowerNumber>();
        implicitly_convertible<FloatValue<PR>,ApproximateNumber>();
        implicitly_convertible<FloatBall<PR>,ValidatedNumber>();
        implicitly_convertible<FloatBall<PR>,ApproximateNumber>();
        implicitly_convertible<FloatBounds<PR>,ValidatedNumber>();
        implicitly_convertible<FloatBounds<PR>,ValidatedUpperNumber>();
        implicitly_convertible<FloatBounds<PR>,ValidatedLowerNumber>();
        implicitly_convertible<FloatBounds<PR>,ApproximateNumber>();
// FIXME: Reintroduce these conversions
//        implicitly_convertible<FloatUpperBound<PR>,ValidatedUpperNumber>();
//        implicitly_convertible<FloatLowerBound<PR>,ValidatedLowerNumber>();
        implicitly_convertible<FloatApproximation<PR>,ApproximateNumber>();
    }
}

/*

namespace Ariadne {
template<class X> struct PythonName<Complex<X>> {
    const char* get() const { static const String res = String("Complex")+class_name<X>(); return res.c_str(); }
};
template<> struct PythonName<Complex<FloatDPApproximation>> {
    const char* get() const { return "ComplexFloatDPApproximation"; }
};
const Complex<Rational> qi = Complex<Rational>(0,1);
} // namespace Ariadne

template<class X> Void export_complex(pymodule& module) {
    typedef Complex<X> Z;
    pybind11::class_<Z> complex_class(module,python_name<Z>());
    complex_class.def(init<X>());
    complex_class.def(init<X,X>());
    if constexpr (HasGenericType<X>) {
        typedef GenericType<X> Y;
        typedef PrecisionType<X> PR;
        complex_class.def(init<Y,Y,PR>());
        complex_class.def(init<Complex<Y>,PR>());
        define_mixed_arithmetic<Complex<X>,Complex<Y>>(module,complex_class);

        // FIXME: Below should not be needed
        complex_class.def(init<Complex<Real>,PR>());
        define_mixed_arithmetic<Complex<X>,Complex<Real>>(module,complex_class);
    }

//    define_mixed_arithmetic<Complex<Integer>,X>(module,complex_class);

    complex_class.def("real_part", &Z::real_part);
    complex_class.def("imaginary_part", &Z::imaginary_part);
    complex_class.def("modulus", &Z::modulus);
    complex_class.def("argument", &Z::argument);

    define_arithmetic<Z>(module,complex_class);
    if constexpr (not Same<X,Rational>) {
        module.def("sqrt", _sqrt_<Z>);
        module.def("exp", _exp_<Z>);
        module.def("log", _log_<Z>);
    }

    module.def("abs", _abs_<Z>);
//    module.def("mag", _mag_<Z>);
    module.def("arg", _arg_<Z>);
    module.def("conj", _conj_<Z>);

    complex_class.def("__str__", &__cstr__<Z>);
//    complex_class.def("__repr__", &__repr__<Z>);

    implicitly_convertible<X,Z>();

    if constexpr (Same<X,Real>) {
        complex_class.def(init<Complex<Integer>>());
        implicitly_convertible<Complex<Integer>,Complex<Real>>();
        complex_class.def(init<Complex<Rational>>());
        implicitly_convertible<Complex<Rational>,Complex<Real>>();
        define_mixed_arithmetic<Complex<Real>,Complex<Integer>>(module,complex_class);
    }
    if constexpr (Same<X,Rational>) {
        complex_class.def(init<Complex<Integer>>());
        implicitly_convertible<Complex<Integer>,Complex<Rational>>();
    }

}

template<> Void export_complex<Integer>(pymodule& module) {
    pybind11::class_<Complex<Integer>> complex_class(module,"ComplexInteger");
    complex_class.def("__str__", &__cstr__<Complex<Integer>>);
    module.attr("i") = Constants::i;
}


Void export_complex_numbers(pymodule& module) {
    export_complex<Integer>(module);
    export_complex<Rational>(module);
    export_complex<Real>(module);
    export_complex<FloatDPBounds>(module);
    export_complex<FloatDPApproximation>(module);

}
*/

/*
namespace Ariadne {
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Value<FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Value<FloatMP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Ball<FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Ball<FloatMP,FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Ball<FloatMP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Bounds<FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Bounds<FloatMP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<UpperBound<FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<UpperBound<FloatMP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<LowerBound<FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<LowerBound<FloatMP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Approximation<FloatDP>> const&);
template OutputStream& operator<<(OutputStream&, PythonRepresentation<Approximation<FloatMP>> const&);
}
*/

Void numeric_submodule(pymodule& module) {
    export_effort(module);
    export_accuracy(module);

    export_logicals(module);

    export_builtins(module);

    export_integer(module);
    export_dyadic(module);
    export_decimal(module);
    export_rational(module);
    export_real(module);

    export_dyadic_bounds(module);
    export_rational_bounds(module);

    export_rounding_mode(module);
    export_precision<DoublePrecision>(module);
    export_precision<MultiplePrecision>(module);
    export_raw_float<DoublePrecision>(module);
    export_raw_float<MultiplePrecision>(module);

    export_user_floats<DoublePrecision>(module);
    export_user_floats<MultiplePrecision>(module);
    export_float_ball<MultiplePrecision,DoublePrecision>(module);

    export_numbers(module);

    export_float_to_number_conversions<DoublePrecision>(module);
    export_float_to_number_conversions<MultiplePrecision>(module);
/*
    export_complex_numbers(module);
*/
}

