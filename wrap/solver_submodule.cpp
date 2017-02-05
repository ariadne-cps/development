/***************************************************************************
 *            solver_submodule.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include "boost_python.hpp"
#include "utilities.hpp"

#include <boost/python.hpp>

#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "solvers/solver_interface.hpp"
#include "solvers/solver.hpp"
#include "function/taylor_function.hpp"

#include "solvers/integrator_interface.hpp"
#include "solvers/integrator.hpp"
#include "solvers/runge_kutta_integrator.hpp"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {

typedef Vector<ValidatedNumericType> ValidatedPointType;
typedef Vector<ApproximateNumericType> ApproximatePointType;

template<class X1, class X2>
Bool operator<(const Vector<X1>& v1, const Vector<X2>& v2);

class SolverWrapper
  : public SolverInterface, public wrapper< SolverInterface >
{
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    Void set_maximum_error(RawFloat64) { this->get_override("set_maximum_error")(); }
    Float64Value maximum_error() const { return this->get_override("maximum_error")(); }
    Void set_maximum_number_of_steps(Nat) { this->get_override("set_maximum_number_of_steps")(); }
    Nat maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    ValidatedPointType zero(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("zero")(); }
    ValidatedPointType fixed_point(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("fixed_point")(); }
    ValidatedPointType solve(const ValidatedVectorFunction& f, const ValidatedPointType& pt) const {
        return this->get_override("solve")(); }
    ValidatedPointType solve(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve")(); }
    ValidatedVectorFunctionModel64 implicit(const ValidatedVectorFunction& f, const ExactBoxType& pd, const ExactBoxType& bx) const {
        return this->get_override("implicit")(); }
    ValidatedScalarFunctionModel64 implicit(const ValidatedScalarFunction& f, const ExactBoxType& pd, const ExactIntervalType& ivl) const {
        return this->get_override("implicit")(); }
    ValidatedVectorFunctionModel64 continuation(const ValidatedVectorFunction& f, const ApproximatePointType& a, const ExactBoxType& X,  const ExactBoxType& A) const {
        return this->get_override("continuation")(); }
    Set< ValidatedPointType > solve_all(const ValidatedVectorFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve_all")(); }
    Void write(OutputStream&) const { this->get_override("write")(); }
};


class IntegratorWrapper
  : public IntegratorInterface, public wrapper< IntegratorInterface >
{
  public:
    IntegratorInterface* clone() const {
        return this->get_override("clone")(); }
    Void set_temporal_order(uint) {
        this->get_override("set_temporal_order")(); }
    Void set_maximum_error(double) {
        this->get_override("set_maximum_error")(); }
    double maximum_error() const {
        return this->get_override("maximum_error")(); }
    Pair<Float64Value,UpperBoxType> flow_bounds(const ValidatedVectorFunction&,const ExactBoxType&,const RawFloat64&) const {
        return this->get_override("flow_bounds")(); }
    ValidatedVectorFunctionModel64 flow_step(const ValidatedVectorFunction&,const ExactBoxType&,RawFloat64&) const {
        return this->get_override("flow_step")(); }
    ValidatedVectorFunctionModel64 flow_step(const ValidatedVectorFunction&,const ExactBoxType&,const Float64Value&,const UpperBoxType&) const {
        return this->get_override("flow_step")(); }
    ValidatedVectorFunctionModel64 flow_to(const ValidatedVectorFunction& vector_field,const ExactBoxType&,const Real&) const {
        return this->get_override("flow_to")(); }
    List<ValidatedVectorFunctionModel64> flow(const ValidatedVectorFunction&,const ExactBoxType&,const Real&,const Real&) const {
        return this->get_override("flow")(); }
    List<ValidatedVectorFunctionModel64> flow(const ValidatedVectorFunction&,const ExactBoxType&,const Real&) const {
        return this->get_override("flow")(); }
    Void write(OutputStream&) const {
        this->get_override("write")(); }
};


} // namespace Ariadne


Void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("solve",pure_virtual((Vector<ValidatedNumericType>(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&)const) &SolverInterface::solve));
    solver_wrapper_class.def("implicit",pure_virtual((ValidatedVectorFunctionModel64(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const ExactBoxType&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("implicit",pure_virtual((ValidatedScalarFunctionModel64(SolverInterface::*)(const ValidatedScalarFunction&,const ExactBoxType&,const ExactIntervalType&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("solve_all",pure_virtual((Set< Vector<ValidatedNumericType> >(SolverInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&)const) &SolverInterface::solve_all));
    //solver_wrapper_class.def(self_ns::str(self));

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());
    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
}



Void export_integrator()
{
    class_<IntegratorWrapper, boost::noncopyable> integrator_wrapper_class("IntegratorInterface");
    integrator_wrapper_class.def("flow_bounds",(Pair<Float64Value,UpperBoxType>(IntegratorInterface::*)(const ValidatedVectorFunction&, const ExactBoxType&, const RawFloat64&)const)&IntegratorInterface::flow_bounds);
    integrator_wrapper_class.def("flow_step",(ValidatedVectorFunctionModel64(IntegratorInterface::*)(const ValidatedVectorFunction&, const ExactBoxType&, RawFloat64&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_step",(ValidatedVectorFunctionModel64(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Float64Value&,const UpperBoxType&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_to",(ValidatedVectorFunctionModel64(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&IntegratorInterface::flow_to);
    integrator_wrapper_class.def("flow",(List<ValidatedVectorFunctionModel64>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&IntegratorInterface::flow);


    class_<TaylorPicardIntegrator > taylor_picard_integrator_class("TaylorPicardIntegrator",init<double>());
    taylor_picard_integrator_class.def("flow_bounds",(Pair<Float64,UpperBoxType>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Float64&)const)&TaylorPicardIntegrator::flow_bounds);
    taylor_picard_integrator_class.def("flow_step", (ValidatedVectorFunctionModel64(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,RawFloat64&)const)&TaylorPicardIntegrator::flow_step);
    taylor_picard_integrator_class.def("flow_step", (ValidatedVectorFunctionModel64(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Float64Value&,const UpperBoxType&)const)&TaylorPicardIntegrator::flow_step);
    taylor_picard_integrator_class.def("flow_to",(ValidatedVectorFunctionModel64(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&TaylorPicardIntegrator::flow_to);
    taylor_picard_integrator_class.def("flow",(List<ValidatedVectorFunctionModel64>(TaylorPicardIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&TaylorPicardIntegrator::flow);

    class_<TaylorSeriesIntegrator > taylor_series_integrator_class("TaylorSeriesIntegrator",init<double>());
    taylor_series_integrator_class.def("maximum_spacial_order",&TaylorSeriesIntegrator::maximum_spacial_order);
    taylor_series_integrator_class.def("maximum_temporal_order",&TaylorSeriesIntegrator::maximum_temporal_order);
    taylor_series_integrator_class.def("maximum_error",&TaylorSeriesIntegrator::maximum_error);
    taylor_series_integrator_class.def("maximum_step_size",&TaylorSeriesIntegrator::maximum_step_size);
    taylor_series_integrator_class.def("set_maximum_spacial_order",&TaylorSeriesIntegrator::set_maximum_spacial_order);
    taylor_series_integrator_class.def("set_maximum_temporal_order",&TaylorSeriesIntegrator::set_maximum_temporal_order);
    taylor_series_integrator_class.def("set_maximum_error",&TaylorSeriesIntegrator::set_maximum_error);
    taylor_series_integrator_class.def("set_maximum_step_size",&TaylorSeriesIntegrator::set_maximum_step_size);
    taylor_series_integrator_class.def("flow_bounds",(Pair<Float64Value,UpperBoxType>(IntegratorInterface::*)(const ValidatedVectorFunction&,const ExactBoxType&,const RawFloat64&)const)&TaylorSeriesIntegrator::flow_bounds);
    taylor_series_integrator_class.def("flow_step", (ValidatedVectorFunctionModel64(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,RawFloat64&)const)&TaylorSeriesIntegrator::flow_step);
    taylor_series_integrator_class.def("flow_step", (ValidatedVectorFunctionModel64(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Float64Value&,const UpperBoxType&)const)&TaylorSeriesIntegrator::flow_step);
    taylor_series_integrator_class.def("flow_to",(ValidatedVectorFunctionModel64(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&TaylorSeriesIntegrator::flow_to);
    taylor_series_integrator_class.def("flow",(List<ValidatedVectorFunctionModel64>(TaylorSeriesIntegrator::*)(const ValidatedVectorFunction&,const ExactBoxType&,const Real&)const)&TaylorSeriesIntegrator::flow);

    class_<RungeKutta4Integrator > runge_kutta_4_integrator_class("RungeKutta4Integrator",init<double>());
    runge_kutta_4_integrator_class.def("step", &RungeKutta4Integrator::step);
    runge_kutta_4_integrator_class.def("evolve", &RungeKutta4Integrator::evolve);
}


Void solver_submodule()
{
    to_python_list< Set< ExactBoxType > >();
    to_python< List< ExactBoxType > >();
    to_python< Pair< Float64, ExactBoxType > >();

    export_solver();
    export_integrator();

}


