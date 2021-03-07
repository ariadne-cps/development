/***************************************************************************
 *            solver_submodule.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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

#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "solvers/solver_interface.hpp"
#include "solvers/solver.hpp"
#include "function/taylor_function.hpp"

#include "solvers/integrator_interface.hpp"
#include "solvers/integrator.hpp"
#include "solvers/runge_kutta_integrator.hpp"

using namespace Ariadne;

namespace Ariadne {

class SolverWrapper
  : public pybind11::wrapper< SolverInterface >
{
    typedef SolverInterface::ValidatedNumericType ValidatedNumericType;
    typedef SolverInterface::ApproximateNumericType ApproximateNumericType;
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    Void set_maximum_error(ApproximateDouble me) { this->get_override("set_maximum_error")(me); }
    ExactDouble maximum_error() const { return this->get_override("maximum_error")(); }
    Void set_maximum_number_of_steps(Nat ns) { this->get_override("set_maximum_number_of_steps")(ns); }
    Nat maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    Vector<ValidatedNumericType> zero(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("zero")(f,bx); }
    Vector<ValidatedNumericType> fixed_point(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("fixed_point")(f,bx); }
    Vector<ValidatedNumericType> solve(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& pt) const {
        return this->get_override("solve")(f,pt); }
    Vector<ValidatedNumericType> solve(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve")(f,bx); }
    ValidatedVectorMultivariateFunctionModelDP implicit(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& pd, const ExactBoxType& bx) const {
        return this->get_override("implicit")(f,pd,bx); }
    ValidatedScalarMultivariateFunctionModelDP implicit(const ValidatedScalarMultivariateFunction& f, const ExactBoxType& pd, const ExactIntervalType& ivl) const {
        return this->get_override("implicit")(f,pd,ivl); }
    ValidatedVectorMultivariateFunctionModelDP continuation(const ValidatedVectorMultivariateFunction& f, const Vector< ApproximateNumericType>& a, const ExactBoxType& X,  const ExactBoxType& A) const {
        return this->get_override("continuation")(f,a,X,A); }
    Set< Vector<ValidatedNumericType> > solve_all(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve_all")(f,bx); }
    Void _write(OutputStream& os) const { this->get_override("_write")(os); }
};


class IntegratorWrapper
  : public pybind11::wrapper<IntegratorInterface>
{
  public:
    IntegratorInterface* clone() const {
        return this->get_override("clone")(); }
    Pair<StepSizeType,UpperBoxType> flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& h) const {
        return this->get_override("flow_bounds")(vf,D,h); }
    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& hst, StepSizeType& hsug) const {
        return this->get_override("flow_step")(vf,D,hst,hsug); }
    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& h, const UpperBoxType& B) const {
        return this->get_override("flow_step")(vf,D,h,B); }
    Void _write(OutputStream& os) const {
        this->get_override("_write")(os); }
};

} // namespace Ariadne


Void export_solvers(pybind11::module& module)
{
    typedef SolverInterface::ValidatedNumericType ValidatedNumericType;
    pybind11::class_<SolverInterface,SolverWrapper> solver_interface_class(module,"SolverInterface");
    solver_interface_class.def("solve", (Vector<ValidatedNumericType>(SolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &SolverInterface::solve);
    solver_interface_class.def("implicit",(ValidatedVectorMultivariateFunctionModelDP(SolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&,const ExactBoxType&)const) &SolverInterface::implicit);
    solver_interface_class.def("implicit",(ValidatedScalarMultivariateFunctionModelDP(SolverInterface::*)(const ValidatedScalarMultivariateFunction&,const ExactBoxType&,const ExactIntervalType&)const) &SolverInterface::implicit);
    solver_interface_class.def("solve_all",(Set< Vector<ValidatedNumericType> >(SolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &SolverInterface::solve_all);
    solver_interface_class.def("__str__",&__cstr__<SolverInterface>);

    pybind11::class_<IntervalNewtonSolver, SolverInterface> interval_newton_solver_class(module,"IntervalNewtonSolver");
    interval_newton_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

    pybind11::class_<KrawczykSolver, SolverInterface> krawczyk_solver_class(module,"KrawczykSolver");
    krawczyk_solver_class.def(pybind11::init<ApproximateDouble,Nat>());
}



Void export_integrators(pybind11::module& module)
{
    pybind11::class_<FlowStepModelType> flow_step_model_class(module,"FlowStepModelType");
    flow_step_model_class.def("__str__", &__cstr__<FlowStepModelType>);
    pybind11::class_<FlowModelType> flow_model_class(module,"FlowModelType");
    flow_model_class.def("size",&FlowModelType::size);
    // Use a lambda here to prevent errors when using Clang/GCC
    flow_model_class.def("__len__",[](FlowModelType const& fm){return fm.size();});
    // NOTE: The export below gives 'ValueError: vector::reserve' at runtime when compiled using Clang
    // flow_model_class.def("__len__",&FlowModelType::size);
    // NOTE: The export below produces 'internal compiler error: in fold_convert_loc' with GCC
    // flow_model_class.def("__len__",(SizeType(FlowModelType::*)()const)&FlowModelType::size);
    flow_model_class.def("__getitem__",&__getitem__<FlowModelType,SizeType>);
    flow_model_class.def("__str__", &__cstr__<FlowModelType>);

    pybind11::class_<IntegratorInterface,IntegratorWrapper> integrator_interface_class(module,"IntegratorInterface");
    integrator_interface_class.def("flow_bounds",(Pair<StepSizeType,UpperBoxType>(IntegratorInterface::*)(const ValidatedVectorMultivariateFunction&, const ExactBoxType&, const StepSizeType&)const)&IntegratorInterface::flow_bounds);
    integrator_interface_class.def("flow_step",(FlowStepModelType(IntegratorInterface::*)(const ValidatedVectorMultivariateFunction&, const ExactBoxType&, const StepSizeType&, StepSizeType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("flow_step",(FlowStepModelType(IntegratorInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&,const StepSizeType&,const UpperBoxType&)const)&IntegratorInterface::flow_step);
    integrator_interface_class.def("__str__", &__cstr__<IntegratorInterface>);
/*
 * TODO: restore
    pybind11::class_<TaylorPicardIntegrator,IntegratorInterface> taylor_picard_integrator_class(module,"TaylorPicardIntegrator");
    taylor_picard_integrator_class.def(pybind11::init<ApproximateDouble>());
    taylor_picard_integrator_class.def("minimum_temporal_order",&TaylorPicardIntegrator::minimum_temporal_order);
    taylor_picard_integrator_class.def("maximum_temporal_order",&TaylorPicardIntegrator::maximum_temporal_order);
    taylor_picard_integrator_class.def("set_minimum_temporal_order",&TaylorPicardIntegrator::set_minimum_temporal_order);
    taylor_picard_integrator_class.def("set_maximum_temporal_order",&TaylorPicardIntegrator::set_maximum_temporal_order);
*/
    pybind11::class_<TaylorSeriesIntegrator,IntegratorInterface> taylor_series_integrator_class(module,"TaylorSeriesIntegrator");
    taylor_series_integrator_class.def(pybind11::init<ApproximateDouble,Nat>());

    pybind11::class_<GradedTaylorSeriesIntegrator,IntegratorInterface> graded_taylor_series_integrator_class(module,"GradedTaylorSeriesIntegrator");
    graded_taylor_series_integrator_class.def(pybind11::init<ApproximateDouble>());

    pybind11::class_<RungeKutta4Integrator> runge_kutta_4_integrator_class(module,"RungeKutta4Integrator");
    runge_kutta_4_integrator_class.def(pybind11::init<ApproximateDouble>());
    runge_kutta_4_integrator_class.def("step", &RungeKutta4Integrator::step);
    runge_kutta_4_integrator_class.def("evolve", &RungeKutta4Integrator::evolve);
}


Void solver_submodule(pybind11::module& module)
{
    export_solvers(module);
    export_integrators(module);
}


