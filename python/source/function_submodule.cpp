/***************************************************************************
 *            function_submodule.cpp
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

#include "utilities.hpp"

#include <iostream>
#include <iomanip>

#include "utility/array.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/differential.hpp"
#include "function/formula.hpp"
#include "function/polynomial.hpp"
#include "function/affine.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/function_expression.hpp"
#include "solvers/constraint_solver.hpp"

#include "function/function_mixin.tpl.hpp"

namespace Ariadne {




MultiIndex multi_index_from_python(pybind11::tuple pytup) {
    MultiIndex res(static_cast<SizeType>(len(pytup)));
    for(SizeType i=0; i!=res.size(); ++i) { res.set(i,pybind11::cast<Nat>(pytup[i])); }
    return res;
}

template<class X, class D>
inline Matrix<X> get_jacobian(const Vector<D>& d) {
    const Nat rs=d.size(); const Nat as=d[0].argument_size();
    Matrix<X> J(rs,as);
    for(Nat i=0; i!=rs; ++i) {
        for(Nat j=0; j!=as; ++j) {
            J[i][j]=d[i][j];
        }
    }
    return J;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<MultiIndex>& arepr) {
    MultiIndex const& a=arepr.reference(); os << "("; for(SizeType i=0; i!=a.size(); ++i) { if(i!=0) { os << ','; } os << a[i]; } os << ")"; return os;
}

template<class P, class SIG> OutputStream& operator<<(OutputStream& os, const Representation< Function<P,SIG> >& frepr) {
    static_cast<const FunctionInterface<P,SIG>&>(frepr.reference()).repr(os); return os;
}

}// namespace Ariadne

using namespace Ariadne;

typedef ExactIntervalType I;
typedef Vector<FloatDPApproximation> FV;
typedef Vector<ExactIntervalType> IV;
typedef Matrix<FloatDPApproximation> FMx;
typedef Matrix<ExactIntervalType> IMx;
typedef Vector< Differential<FloatDPApproximation> > FSDV;
typedef Vector< Differential<ExactIntervalType> > ISDV;
typedef Vector<ValidatedTaylorModelDP> TMV;
typedef ValidatedVectorMultivariateTaylorFunctionModelDP TFM;
typedef ValidatedTaylorModelDP TM;


template<class AT, class F> void define_call(pybind11::class_<F>& f) {
    typedef typename F::template Argument<AT> A; typedef ResultOf<F(A)> R; f.def("__call__", (R(F::*)(A const&)const) &F::operator()); }

template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<ScalarUnivariateFunction<P>>& function_class) {
    typedef ScalarUnivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(slope(declval<F>(),declval<A>())) R;
    module.def("slope", (R(*)(F const&,A const&)) &slope);
    function_class.def("slope", (R(F::*)(A const&)const) &F::slope);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<ScalarMultivariateFunction<P>>& function_class) {
    typedef ScalarMultivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(gradient(declval<F>(),declval<A>())) R;
    module.def("gradient", (R(*)(F const&,A const&)) &gradient);
    function_class.def("gradient", (R(F::*)(A const&)const) &F::gradient);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<VectorUnivariateFunction<P>>& function_class) {
    typedef VectorUnivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(tangent(declval<F>(),declval<A>())) R;
    module.def("tangent", (R(*)(F const&,A const&)) &tangent);
    function_class.def("tangent", (R(F::*)(A const&)const) &F::tangent);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<VectorMultivariateFunction<P>>& function_class) {
    typedef VectorMultivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(jacobian(declval<F>(),declval<A>())) R;
    module.def("jacobian", (R(*)(F const&,A const&)) &jacobian);
    function_class.def("jacobian", (R(F::*)(A const&)const) &F::jacobian);
}


template<class F, class T> using ArgumentType = typename F::template Argument<T>;

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, ApproximateTag) {
    define_call<FloatDPApproximation>(function_class);
    define_call<FloatMPApproximation>(function_class);
    define_call<Differential<FloatDPApproximation>>(function_class);
    define_call<Differential<FloatMPApproximation>>(function_class);

    function_class.def("differential", _differential_<F,ArgumentType<F,FloatDPApproximation>,DegreeType>);
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatMPApproximation>,DegreeType>);

    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatDPApproximation>>);
    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatMPApproximation>>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatDPApproximation>,DegreeType>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatMPApproximation>,DegreeType>);

    define_named_derivatives<FloatDPApproximation>(module,function_class);
    define_named_derivatives<FloatMPApproximation>(module,function_class);
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, ValidatedTag) {
    define_call<FloatDPBounds>(function_class);
    define_call<FloatMPBounds>(function_class);
    define_call<Differential<FloatDPBounds>>(function_class);
    define_call<Differential<FloatMPBounds>>(function_class);

    function_class.def("differential", _differential_<F,ArgumentType<F,FloatDPBounds>,DegreeType>);
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatMPBounds>,DegreeType>);

    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatDPBounds>>);
    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatMPBounds>>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatDPBounds>,DegreeType>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatMPBounds>,DegreeType>);

    define_named_derivatives<FloatDPBounds>(module,function_class);
    define_named_derivatives<FloatMPBounds>(module,function_class);

    export_function_evaluation(module,function_class,ApproximateTag());
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, EffectiveTag) {
    export_function_evaluation(module,function_class,ValidatedTag());
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class)
{
    export_function_evaluation(module, function_class, Paradigm<F>());
}


Void export_multi_index(pybind11::module& module)
{
    pybind11::class_< MultiIndex > multi_index_class(module,"MultiIndex");
    multi_index_class.def(pybind11::init<Nat>());
    multi_index_class.def(pybind11::init(&multi_index_from_python));
    multi_index_class.def(pybind11::init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def("__str__", &__cstr__<MultiIndex>);
    multi_index_class.def("__repr__", &__repr__<MultiIndex>);

    pybind11::implicitly_convertible<pybind11::tuple,MultiIndex>();
}


template<class X>
Void export_polynomial(pybind11::module& module)
{
    pybind11::class_< MultivariateMonomial<X> > monomial_class(module,python_name<X>("MultivariateMonomial").c_str());
    monomial_class.def(pybind11::init<MultiIndex,X>());
    monomial_class.def("index", (MultiIndex const&(MultivariateMonomial<X>::*)()const) &MultivariateMonomial<X>::index);
    monomial_class.def("coefficient", (X const&(MultivariateMonomial<X>::*)()const) &MultivariateMonomial<X>::coefficient);
    monomial_class.def("__str__", &__cstr__<MultivariateMonomial<X>>);


    pybind11::class_< MultivariatePolynomial<X> > polynomial_class(module,python_name<X>("MultivariatePolynomial").c_str());
    polynomial_class.def(pybind11::init< MultivariatePolynomial<X> >());
    polynomial_class.def_static("constant", (MultivariatePolynomial<X>(*)(SizeType,X const&)) &MultivariatePolynomial<X>::constant);

    if constexpr (HasPrecisionType<X>) {
        typedef typename X::PrecisionType PR;
        polynomial_class.def(pybind11::init<Nat,PR>());
        polynomial_class.def_static("variable", (MultivariatePolynomial<X>(*)(SizeType,SizeType,PR)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("coordinate", (MultivariatePolynomial<X>(*)(SizeType,SizeType,PR)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("variables", [](Nat as,PR pr){return MultivariatePolynomial<X>::variables(as,pr).array();});
    } else {
        polynomial_class.def(pybind11::init<Nat>());
        polynomial_class.def_static("variable", (MultivariatePolynomial<X>(*)(SizeType,SizeType)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("coordinate", (MultivariatePolynomial<X>(*)(SizeType,SizeType)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("variables", [](Nat as){return MultivariatePolynomial<X>::variables(as).array();});
    }

    polynomial_class.def("argument_size", &MultivariatePolynomial<X>::argument_size);
    polynomial_class.def("insert", &MultivariatePolynomial<X>::insert);

    define_algebra(module,polynomial_class);
    polynomial_class.def("__str__",&__cstr__<MultivariatePolynomial<X>>);

    export_vector<MultivariatePolynomial<X>>(module, (python_name<X>("MultivariatePolynomialVector")).c_str());
}

Void export_polynomials(pybind11::module& module)
{
    export_polynomial<FloatDPBounds>(module);
    export_polynomial<FloatDPApproximation>(module);
}

Void export_domains(pybind11::module& module)
{
    pybind11::class_<RealDomain> real_domain_class(module,"RealDomain");
    real_domain_class.def(pybind11::init<>());
    //    real_domain_class.def(pybind11::init<SizeOne>());
    real_domain_class.def("__repr__", &__cstr__<RealDomain>);

    pybind11::class_<EuclideanDomain> euclidean_domain_class(module,"EuclideanDomain");
    euclidean_domain_class.def(pybind11::init<SizeType>());
    euclidean_domain_class.def("__repr__", &__cstr__<EuclideanDomain>);
}

template<class P> Void export_scalar_univariate_function(pybind11::module& module)
{
    pybind11::class_<ScalarUnivariateFunction<P>> function_class(module,(class_name<P>()+"ScalarUnivariateFunction").c_str());
    function_class.def(pybind11::init<ScalarUnivariateFunction<P>>());
    if constexpr (Same<P,ValidatedTag>) {
        function_class.def(pybind11::init<ScalarUnivariateFunction<EffectiveTag>>());
//        pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ValidatedTag>>();
    }
    function_class.def("derivative", (ScalarUnivariateFunction<P>(ScalarUnivariateFunction<P>::*)(SizeOne)const) &ScalarUnivariateFunction<P>::derivative);
    function_class.def("derivative", (ScalarUnivariateFunction<P>(*)(const ScalarUnivariateFunction<P>&)) &derivative);

    export_function_evaluation(module,function_class);
    if constexpr (not Same<P,ApproximateTag>) {
        module.def("derivative", [](const ScalarUnivariateFunction<P>& f, const FloatDPBounds& x){return static_cast<FloatDPBounds>(differential(f,x,1u).gradient()[0]);} );
        module.def("second_derivative", [](const ScalarUnivariateFunction<P>& f, const FloatDPBounds& x){return static_cast<FloatDPBounds>(differential(f,x,2u).hessian()[0][0]);} );
    }
    module.def("derivative", [](const ScalarUnivariateFunction<P>& f, const FloatDPApproximation& x){return static_cast<FloatDPApproximation>(differential(f,x,1u).gradient()[0]);} );
    module.def("second_derivative", [](const ScalarUnivariateFunction<P>& f, const FloatDPApproximation& x){return static_cast<FloatDPApproximation>(differential(f,x,2u).hessian()[0][0]);} );

    function_class.def_static("constant", (ScalarUnivariateFunction<P>(*)(RealDomain,Number<P>)) &ScalarUnivariateFunction<P>::constant);
    function_class.def_static("coordinate", (ScalarUnivariateFunction<P>(*)()) &ScalarUnivariateFunction<P>::coordinate);
    function_class.def_static("coordinate", (ScalarUnivariateFunction<P>(*)(RealDomain)) &ScalarUnivariateFunction<P>::coordinate);
    function_class.def_static("identity", (ScalarUnivariateFunction<P>(*)()) &ScalarUnivariateFunction<P>::identity);
    function_class.def_static("identity", (ScalarUnivariateFunction<P>(*)(RealDomain)) &ScalarUnivariateFunction<P>::identity);

    define_elementary_algebra<ScalarUnivariateFunction<P>,Number<P>>(module,function_class);

    function_class.def("__str__", &__cstr__<ScalarUnivariateFunction<P>>);
}

template<class P> Void export_vector_univariate_function(pybind11::module& module)
{
    pybind11::class_<VectorUnivariateFunction<P>> vector_univariate_function_class(module,(class_name<P>()+"VectorUnivariateFunction").c_str());
    vector_univariate_function_class.def(pybind11::init<VectorUnivariateFunction<P>>());
    vector_univariate_function_class.def(pybind11::init<SizeType,RealDomain>());
    if constexpr (Same<P,ValidatedTag>) {
        vector_univariate_function_class.def(pybind11::init<VectorUnivariateFunction<EffectiveTag>>());
//        pybind11::implicitly_convertible<VectorUnivariateFunction<EffectiveTag>,VectorUnivariateFunction<ValidatedTag>>();
    }
    // NOTE: This must go *after* the conversion constructor
    vector_univariate_function_class.def(pybind11::init([](std::vector<ScalarUnivariateFunction<P>> const& lst){return VectorUnivariateFunction<P>(lst);}));

    vector_univariate_function_class.def("result_size", &VectorUnivariateFunction<P>::result_size);
    vector_univariate_function_class.def("argument_size", &VectorUnivariateFunction<P>::argument_size);
    vector_univariate_function_class.def("__getitem__", &VectorUnivariateFunction<P>::get);
    vector_univariate_function_class.def("__setitem__", &VectorUnivariateFunction<P>::set);

    // TODO: Put these in C++ API
    // define_vector_algebra_arithmetic<VectorUnivariateFunction<P>,ScalarUnivariateFunction<P>,Number<P>>(module,vector_univariate_function_class);

    // FIXME: Define vector function operations for Validated and Approximate
//    if constexpr (Same<P,EffectiveTag>) {
//        define_vector_algebra_arithmetic<VectorUnivariateFunction<P>,ScalarUnivariateFunction<P>>(module,vector_univariate_function_class);
//    }
    export_function_evaluation(module,vector_univariate_function_class);

    vector_univariate_function_class.def("__str__", &__cstr__<VectorUnivariateFunction<P>>);
    vector_univariate_function_class.def("__repr__", &__crepr__<VectorUnivariateFunction<P>>);

//    module.def("join", (VectorUnivariateFunction<P>(*)(const ScalarUnivariateFunction<P>&, const ScalarUnivariateFunction<P>&)) &join);
//    module.def("join", (VectorUnivariateFunction<P>(*)(const VectorUnivariateFunction<P>&, const ScalarUnivariateFunction<P>&)) &join);
//    module.def("join", (VectorUnivariateFunction<P>(*)(const ScalarUnivariateFunction<P>&, const VectorUnivariateFunction<P>&)) &join);
//    module.def("join", (VectorUnivariateFunction<P>(*)(const VectorUnivariateFunction<P>&, const VectorUnivariateFunction<P>&)) &join);

    module.def("compose", (ScalarUnivariateFunction<P>(*)(const ScalarUnivariateFunction<P>&,const ScalarUnivariateFunction<P>&)) &compose);
    module.def("compose", (ScalarUnivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&,const VectorUnivariateFunction<P>&)) &compose);
    module.def("compose", (VectorUnivariateFunction<P>(*)(const VectorUnivariateFunction<P>&,const ScalarUnivariateFunction<P>&)) &compose);
    module.def("compose", (VectorUnivariateFunction<P>(*)(const VectorMultivariateFunction<P>&,const VectorUnivariateFunction<P>&)) &compose);
}


Void export_scalar_univariate_functions(pybind11::module& module)
{
    export_scalar_univariate_function<EffectiveTag>(module);
    export_scalar_univariate_function<ValidatedTag>(module);
    export_scalar_univariate_function<ApproximateTag>(module);
    pybind11::implicitly_convertible<ScalarUnivariateFunction<EffectiveTag>,ScalarUnivariateFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<ScalarUnivariateFunction<EffectiveTag>,ScalarUnivariateFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<ScalarUnivariateFunction<ValidatedTag>,ScalarUnivariateFunction<ApproximateTag>>();
}

Void export_vector_univariate_functions(pybind11::module& module)
{
    export_vector_univariate_function<EffectiveTag>(module);
    export_vector_univariate_function<ValidatedTag>(module);
    export_vector_univariate_function<ApproximateTag>(module);
    pybind11::implicitly_convertible<VectorUnivariateFunction<EffectiveTag>,VectorUnivariateFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<VectorUnivariateFunction<EffectiveTag>,VectorUnivariateFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<VectorUnivariateFunction<ValidatedTag>,VectorUnivariateFunction<ApproximateTag>>();
}



template<class P> Void export_scalar_function(pybind11::module& module)
{
    pybind11::class_<ScalarMultivariateFunction<P>> scalar_function_class(module,(class_name<P>()+"ScalarMultivariateFunction").c_str());
    scalar_function_class.def(pybind11::init<ScalarMultivariateFunction<P>>());
    scalar_function_class.def(pybind11::init<SizeType>());
    scalar_function_class.def("argument_size", &ScalarMultivariateFunction<P>::argument_size);
    scalar_function_class.def("derivative", &ScalarMultivariateFunction<P>::derivative);
    define_elementary_algebra<ScalarMultivariateFunction<P>,Number<P>>(module,scalar_function_class);

    if constexpr (Same<P,ValidatedTag>) {
        scalar_function_class.def(pybind11::init<ScalarMultivariateFunction<EffectiveTag>>());
        pybind11::implicitly_convertible<ScalarMultivariateFunction<EffectiveTag>,ScalarMultivariateFunction<ValidatedTag>>();
    }

    if constexpr (Same<P,EffectiveTag>) {
        scalar_function_class.def("__eq__", &__eq__<EffectiveScalarMultivariateFunction,EffectiveNumber,Return<EffectiveConstraint>>);
        scalar_function_class.def("__le__", &__le__<EffectiveScalarMultivariateFunction,EffectiveNumber,Return<EffectiveConstraint>>);
        scalar_function_class.def("__ge__", &__ge__<EffectiveScalarMultivariateFunction,EffectiveNumber,Return<EffectiveConstraint>>);
    }

    scalar_function_class.def("__str__", &__cstr__<ScalarMultivariateFunction<P>>);
    scalar_function_class.def("__repr__", &__crepr__<ScalarMultivariateFunction<P>>);

    scalar_function_class.def_static("constant", (ScalarMultivariateFunction<P>(*)(SizeType,Number<P>)) &ScalarMultivariateFunction<P>::constant);
    scalar_function_class.def_static("coordinate", (ScalarMultivariateFunction<P>(*)(EuclideanDomain,SizeType)) &ScalarMultivariateFunction<P>::coordinate);
    scalar_function_class.def_static("coordinate", (ScalarMultivariateFunction<P>(*)(SizeType,SizeType)) &ScalarMultivariateFunction<P>::coordinate);

    scalar_function_class.def("gradient", (Covector<FloatDPApproximation>(ScalarMultivariateFunction<P>::*)(const Vector<FloatDPApproximation>&)const) &ScalarMultivariateFunction<P>::gradient);
    if constexpr (not Same<P,ApproximateTag>) {
        scalar_function_class.def("gradient", (Covector<FloatDPBounds>(ScalarMultivariateFunction<P>::*)(const Vector<FloatDPBounds>&)const) &ScalarMultivariateFunction<P>::gradient);
    }

    module.def("derivative", (ScalarMultivariateFunction<P>(ScalarMultivariateFunction<P>::*)(SizeType)const) &ScalarMultivariateFunction<P>::derivative);

    module.def("evaluate", (Scalar<FloatDPApproximation>(*)(const ScalarMultivariateFunction<P>&,const Vector<FloatDPApproximation>&)) &evaluate);
    if constexpr (not Same<P,ApproximateTag>) {
        module.def("evaluate", (Scalar<FloatDPBounds>(*)(const ScalarMultivariateFunction<P>&,const Vector<FloatDPBounds>&)) &evaluate);
    }

    export_function_evaluation(module,scalar_function_class);
}


template<class P> Void export_vector_function(pybind11::module& module)
{
    pybind11::class_<VectorMultivariateFunction<P>> vector_function_class(module,(class_name<P>()+"VectorMultivariateFunction").c_str());
    vector_function_class.def(pybind11::init<VectorMultivariateFunction<P>>());
    vector_function_class.def(pybind11::init<SizeType,SizeType>());
    if constexpr (Same<P,ValidatedTag>) {
        vector_function_class.def(pybind11::init<VectorMultivariateFunction<EffectiveTag>>());
//        pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ValidatedTag>>();
    }
    // NOTE: This must go *after* the conversion constructor
    vector_function_class.def(pybind11::init([](std::vector<ScalarMultivariateFunction<P>> const& lst){return VectorMultivariateFunction<P>(lst);}));

    vector_function_class.def("result_size", &VectorMultivariateFunction<P>::result_size);
    vector_function_class.def("argument_size", &VectorMultivariateFunction<P>::argument_size);
    vector_function_class.def("__getitem__", &VectorMultivariateFunction<P>::get);
    vector_function_class.def("__setitem__", &VectorMultivariateFunction<P>::set);

    // TODO: Put these in C++ API
    // define_vector_algebra_arithmetic<VectorMultivariateFunction<P>,ScalarMultivariateFunction<P>,Number<P>>(module,vector_function_class);

    // FIXME: Define vector function operations for Validated and Approximate
    if constexpr (Same<P,EffectiveTag>) {
        define_vector_algebra_arithmetic<VectorMultivariateFunction<P>,ScalarMultivariateFunction<P>>(module,vector_function_class);
    }
    export_function_evaluation(module,vector_function_class);

    vector_function_class.def("__str__", &__cstr__<VectorMultivariateFunction<P>>);
    vector_function_class.def("__repr__", &__crepr__<VectorMultivariateFunction<P>>);

    vector_function_class.def_static("zeros", (VectorMultivariateFunction<P>(*)(SizeType,SizeType)) &VectorMultivariateFunction<P>::zeros);
    vector_function_class.def_static("zeros", (VectorMultivariateFunction<P>(*)(SizeType,EuclideanDomain)) &VectorMultivariateFunction<P>::zeros);
    vector_function_class.def_static("identity", (VectorMultivariateFunction<P>(*)(SizeType)) &VectorMultivariateFunction<P>::identity);
    vector_function_class.def_static("identity", (VectorMultivariateFunction<P>(*)(EuclideanDomain)) &VectorMultivariateFunction<P>::identity);

    module.def("join", (VectorMultivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&, const ScalarMultivariateFunction<P>&)) &join);
    module.def("join", (VectorMultivariateFunction<P>(*)(const VectorMultivariateFunction<P>&, const ScalarMultivariateFunction<P>&)) &join);
    module.def("join", (VectorMultivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&, const VectorMultivariateFunction<P>&)) &join);
    module.def("join", (VectorMultivariateFunction<P>(*)(const VectorMultivariateFunction<P>&, const VectorMultivariateFunction<P>&)) &join);

    module.def("compose", (ScalarMultivariateFunction<P>(*)(const ScalarUnivariateFunction<P>&,const ScalarMultivariateFunction<P>&)) &compose);
    module.def("compose", (ScalarMultivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&,const VectorMultivariateFunction<P>&)) &compose);
    module.def("compose", (VectorMultivariateFunction<P>(*)(const VectorUnivariateFunction<P>&,const ScalarMultivariateFunction<P>&)) &compose);
    module.def("compose", (VectorMultivariateFunction<P>(*)(const VectorMultivariateFunction<P>&,const VectorMultivariateFunction<P>&)) &compose);
}

template<class Y, class X> Void export_procedure(pybind11::module& module) {
    typedef Paradigm<Y> P;
    pybind11::class_<Procedure<Y>> procedure_class(module,(class_name<P>()+"Procedure").c_str());
    procedure_class.def("__str__", &__cstr__<Procedure<Y>>);
    module.def("make_procedure", (Procedure<Y>(*)(ScalarMultivariateFunction<P> const&)) &make_procedure);
    module.def("evaluate", (X(*)(Procedure<Y> const&, Vector<X> const&)) &evaluate);
    module.def("gradient", (Covector<X>(*)(Procedure<Y> const&, Vector<X> const&)) &gradient);
    module.def("hessian", (X(*)(Procedure<Y> const&, Vector<X> const&, Vector<X> const&)) &hessian);
}

Void export_scalar_functions(pybind11::module& module) {
    export_scalar_function<ApproximateTag>(module);
    export_scalar_function<ValidatedTag>(module);
    export_scalar_function<EffectiveTag>(module);
    pybind11::implicitly_convertible<ScalarMultivariateFunction<EffectiveTag>,ScalarMultivariateFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<ScalarMultivariateFunction<EffectiveTag>,ScalarMultivariateFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<ScalarMultivariateFunction<ValidatedTag>,ScalarMultivariateFunction<ApproximateTag>>();
    module.def("lie_derivative", (ScalarMultivariateFunction<EffectiveTag>(*)(const ScalarMultivariateFunction<EffectiveTag>&,const VectorMultivariateFunction<EffectiveTag>&)) &lie_derivative);
}

Void export_vector_functions(pybind11::module& module) {
    export_vector_function<ApproximateTag>(module);
    export_vector_function<ValidatedTag>(module);
    export_vector_function<EffectiveTag>(module);
    pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<VectorMultivariateFunction<ValidatedTag>,VectorMultivariateFunction<ApproximateTag>>();
}




Void function_submodule(pybind11::module& module) {

    export_multi_index(module);

    export_polynomials(module);

    export_domains(module);

    export_scalar_univariate_functions(module);
    export_vector_univariate_functions(module);
    export_scalar_functions(module);
    export_vector_functions(module);

    export_procedure<ApproximateNumber, FloatDPApproximation>(module);
    export_procedure<ValidatedNumber, FloatDPBounds>(module);
}

