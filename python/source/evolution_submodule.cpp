/***************************************************************************
 *            evolution_submodule.cpp
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

#include "dynamics/orbit.hpp"
#include "dynamics/map_evolver.hpp"
#include "dynamics/vector_field_evolver.hpp"

using namespace Ariadne;


template<class ORB>
Void export_orbit(pybind11::module& module, const char* name)
{
    pybind11::class_<ORB> orbit_class(module,name);
    orbit_class.def("reach", &ORB::reach);
    orbit_class.def("evolve", &ORB::final);
    orbit_class.def("final", &ORB::final);
    orbit_class.def("__str__", &__cstr__<ORB>);
}


Void export_semantics(pybind11::module& module) {
    pybind11::enum_<Semantics> semantics_enum(module,"Semantics");
    semantics_enum.value("UPPER", Semantics::UPPER);
    semantics_enum.value("LOWER", Semantics::LOWER);
    //    semantics_enum.def("__repr__" , &__cstr__<Semantics>);
}


template<class SIM> Void export_simulator(pybind11::module& module, const char* name);


template<class EV>
Void export_evolver_interface(pybind11::module& module, const char* name)
{
    pybind11::class_<EV> evolver_interface_class(module,name);
}

template<class EV, class... Params>
Void export_evolver(pybind11::module& module, const char* name)
{
    typedef EV Evolver;
    typedef typename EV::Interface Interface;
    typedef typename EV::EnclosureType EnclosureType;
    typedef typename EV::TerminationType TerminationType;
    typedef typename EV::OrbitType OrbitType;

    pybind11::class_<Evolver,pybind11::bases<Interface>> evolver_class(module,name);
    evolver_class.def(pybind11::init<Params...>());
    evolver_class.def("orbit",(OrbitType(Evolver::*)(const EnclosureType&,const TerminationType&,Semantics)const) &Evolver::orbit);
    evolver_class.def("__str__",&__cstr__<Evolver>);
}


template<class RA, class... Params>
Void export_reachability_analyser(pybind11::module& module, const char* name)
{
    typedef typename RA::ConfigurationType Configuration;
    typedef typename RA::StorageType StorageType;
    typedef typename RA::OvertSetInterfaceType OvertSetType;
    typedef typename RA::CompactSetInterfaceType CompactSetType;
    typedef typename RA::TimeType TimeType;

    auto const& reference_internal = pybind11::return_value_policy::reference_internal;

    pybind11::class_<RA> reachability_analyser_class(module,name);
    reachability_analyser_class.def(pybind11::init<Params...>());
    reachability_analyser_class.def("configuration",(Configuration&(RA::*)())&RA::configuration,reference_internal);
    reachability_analyser_class.def("evolver",&RA::evolver,reference_internal);
    reachability_analyser_class.def("lower_reach",(StorageType(RA::*)(OvertSetType const&,TimeType const&)const) &RA::lower_reach);
    reachability_analyser_class.def("upper_reach", (StorageType(RA::*)(CompactSetType const&,TimeType const&)const) &RA::upper_reach);
    reachability_analyser_class.def("outer_chain_reach", (StorageType(RA::*)(CompactSetType const&)const)&RA::outer_chain_reach);
    reachability_analyser_class.def("__str__",&__cstr__<RA>);
}


Void evolution_submodule(pybind11::module& module)
{
    export_semantics(module);

    export_evolver_interface<MapEvolver::Interface>(module,"MapEvolverInterface");
    export_evolver_interface<VectorFieldEvolver::Interface>(module,"VectorFieldEvolverInterface");

    export_evolver<MapEvolver, IteratedMap>(module,"MapEvolver");
    export_evolver<VectorFieldEvolver, VectorField const&, Configuration<VectorFieldEvolver> const&>(module,"VectorFieldEvolver");
}
