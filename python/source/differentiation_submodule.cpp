/***************************************************************************
 *            differentiation_submodule.cpp
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
#include <pybind11/stl.h>

#include "utilities.hpp"

#include "utility/typedefs.hpp"
#include "utility/array.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/differential.hpp"
#include "algebra/univariate_differential.hpp"
#include "algebra/fixed_differential.hpp"
#include "algebra/fixed_univariate_differential.hpp"
#include "algebra/expansion.inl.hpp"

using namespace Ariadne;

namespace Ariadne {

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Expansion<MultiIndex,X>>& repr);

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Differential<X> >& repr) {
    const Differential<X>& diff=repr.reference();
    os << python_name<X>("Differential").c_str() << "(" << python_representation(diff.expansion()) << "," << diff.degree() << ")";
    return os;
}

} // namespace Ariadne


template<class D> using ValueType = decltype(declval<D>().value());
template<class D> using GradientType = decltype(declval<D>().gradient());
template<class D> using HessianType = decltype(declval<D>().hessian());

template<class DIFF>
Void export_differential(pybind11::module& module, const String& name)
{
    typedef typename DIFF::ValueType X;
    typedef DIFF D;

    pybind11::class_<D> differential_class(module,name.c_str());
    differential_class.def(pybind11::init<D>());

    if constexpr (HasGenericType<X>) {
        typedef typename X::PrecisionType PR;
        differential_class.def( pybind11::init<SizeType,DegreeType,PR>());
    } else {
        differential_class.def( pybind11::init<SizeType,DegreeType>());
    }

    differential_class.def_static("constant",(D(*)(SizeType, DegreeType, const X&)) &D::constant);
    differential_class.def_static("variable",(D(*)(SizeType as, DegreeType deg, const X& y, SizeType j)) &D::variable);
    differential_class.def_static("constants",(Vector<D>(*)(SizeType, DegreeType, const Vector<X>&)) &D::constants);
    differential_class.def_static("variables",(Vector<D>(*)(DegreeType, const Vector<X>&)) &D::variables);
    if constexpr (HasGenericType<X>) {
        typedef typename X::GenericType Y; typedef typename X::PrecisionType PR;
        differential_class.def_static("constant",(D(*)(SizeType, DegreeType, const Y&, const PR&)) &D::constant);
        differential_class.def_static("variable",(D(*)(SizeType as, DegreeType deg, const Y& y, SizeType j, const PR& pr)) &D::variable);
        differential_class.def_static("constants",(Vector<D>(*)(SizeType, DegreeType, const Vector<Y>&, const PR&)) &D::constants);
        differential_class.def_static("variables",(Vector<D>(*)(DegreeType, const Vector<Y>&, const PR&)) &D::variables);
    }

    differential_class.def( pybind11::init<Expansion<MultiIndex,X>,DegreeType>());
    differential_class.def("__getitem__", &__getitem__<D,MultiIndex,X>);
    differential_class.def("__setitem__",&__setitem__<D,MultiIndex,X>);

    define_elementary_algebra<D,X>(module, differential_class);
    define_inplace_algebra<D,X>(module,differential_class);
    differential_class.def("__str__", &__cstr__<D>);
    differential_class.def("__repr__", &__repr__<D>);

    differential_class.def("value", (ValueType<D>(D::*)()const)&D::value);
    differential_class.def("gradient", (GradientType<D>(D::*)()const)&D::gradient);
    differential_class.def("hessian", (HessianType<D>(D::*)()const)&D::hessian);
    differential_class.def("expansion", (Expansion<MultiIndex,X>const&(D::*)()const)&D::expansion);

    differential_class.def_static("constant",(D(*)(SizeType, DegreeType, const X&))&D::constant);
    differential_class.def_static("variable",(D(*)(SizeType, DegreeType, const X&, SizeType))&D::variable);
    differential_class.def_static("constants",(Vector<D>(*)(SizeType, DegreeType, const Vector<X>&))&D::constants);
    differential_class.def_static("variables",(Vector<D>(*)(DegreeType, const Vector<X>&))&D::variables);

    if constexpr (HasGenericType<X>) {
        typedef typename X::GenericType Y;
        define_mixed_arithmetic<D,Y>(module, differential_class);
    }

    module.def("derivative", (D(*)(const D&, SizeType))&D::_derivative);
    module.def("antiderivative", (D(*)(const D&, SizeType))&D::_antiderivative);

}

template<class DIFF>
Void
export_differential_vector(pybind11::module& module, const String& name)
{
    typedef typename DIFF::ValueType X;
    typedef Vector<X> V;
    typedef DIFF D;
    typedef Vector<D> DV;

    pybind11::class_<DV> differential_vector_class(module,name.c_str());
    differential_vector_class.def(pybind11::init<DV>());
    differential_vector_class.def(pybind11::init<Vector<D>>());
    if constexpr (HasPrecisionType<X>) {
        typedef typename X::PrecisionType PR;
        differential_vector_class.def(pybind11::init<Nat,Nat,Nat,PR>());
    } else {
        differential_vector_class.def(pybind11::init<Nat,Nat,Nat>());
    }
    differential_vector_class.def("__getitem__", &__getitem2__<DV,Nat,MultiIndex,X>);
    differential_vector_class.def("__getitem__", &__getitem__<DV,Nat,D>);
    differential_vector_class.def("__setitem__",&__setitem__<DV,Nat,X>);
    differential_vector_class.def("__setitem__",&__setitem__<DV,Nat,D>);

    define_vector_algebra_arithmetic(module, differential_vector_class);

    differential_vector_class.def("size", &DV::size);
    differential_vector_class.def("argument_size", &DV::argument_size);
    differential_vector_class.def("degree", &DV::degree);
    differential_vector_class.def("value", &DV::value);
    differential_vector_class.def("jacobian", &DV::jacobian);
    differential_vector_class.def("__str__",&__cstr__<DV>);
    //differential_vector_class.def("__repr__",&__repr__<DV>);

    module.def("compose",(D(*)(const D&,const DV&))&D::_compose);
    module.def("compose",(DV(*)(const DV&,const DV&))&DV::_compose);

    module.def("solve",(DV(*)(const DV&,const V&))&DV::_solve);
    module.def("flow",(DV(*)(const DV&,const V&))&DV::_flow);

    //NOTE: Not in C++ interface
    //module.def("lie_derivative", (DV(*)(const DV&,const DV&))&lie_derivative);
}

Void differentiation_submodule(pybind11::module& module)
{

    export_differential< Differential<FloatDPApproximation> >(module,python_name<FloatDPApproximation>("Differential"));
    export_differential< Differential<FloatDPBounds> >(module,python_name<FloatDPBounds>("Differential"));
    export_differential_vector< Differential<FloatDPApproximation> >(module,python_name<FloatDPApproximation>("DifferentialVector"));
    export_differential_vector< Differential<FloatDPBounds> >(module,python_name<FloatDPBounds>("DifferentialVector"));

    export_differential< Differential<FloatMPApproximation> >(module,python_name<FloatMPApproximation>("Differential"));
    export_differential< Differential<FloatMPBounds> >(module,python_name<FloatMPBounds>("Differential"));
    export_differential_vector< Differential<FloatMPApproximation> >(module,python_name<FloatMPApproximation>("DifferentialVector"));
    export_differential_vector< Differential<FloatMPBounds> >(module,python_name<FloatMPBounds>("DifferentialVector"));
}

