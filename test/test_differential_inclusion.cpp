/***************************************************************************
 *            test_differential_inclusion.cpp
 *
 *  Copyright  2008-17  Pieter Collins, Sanja Zivanovic
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

#include "dynamics/differential_inclusion.hpp"

#include "algebra/sweeper.hpp"
#include "solvers/integrator_interface.hpp"
#include "geometry/box.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "geometry/function_set.hpp"
#include "output/graphics.hpp"

#include "test/test.hpp"
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

namespace Ariadne {


template<class F, class S> List<ResultOf<F(S)>> map(F const& f, List<S> const& list) {
    List<ResultOf<F(S)>> result; for(auto item : list) { result.append(f(item)); } return result;
}

ValidatedConstrainedImageSet range(ValidatedVectorFunctionModelType const& fm) {
    return ValidatedConstrainedImageSet(fm.domain(),fm);
}

ThresholdSweeperDP make_threshold_sweeper(double thr) { return ThresholdSweeperDP(DoublePrecision(),thr); }
GradedSweeperDP make_graded_sweeper(SizeType deg) { return GradedSweeperDP(DoublePrecision(),deg); }
GradedThresholdSweeperDP make_graded_threshold_sweeper(SizeType deg, double thr) { return GradedThresholdSweeperDP(DoublePrecision(),deg, thr); }

template<class C> struct Reverse {
    C const& _c;
    Reverse(C const& c) :  _c(c) {}
    typename C::const_reverse_iterator begin() const{ return _c.rbegin(); }
    typename C::const_reverse_iterator end() const { return _c.rend(); }
};
template<class C> Reverse<C> reverse(C const& c) { return Reverse<C>(c); }

} // namespace Ariadne

using namespace Ariadne;

class TestInclusionIntegrator {

    Void run_battery_fixed_avg(String name,
                                ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                                RealBox real_starting_set, Real evolution_time, double step, SizeType avg, SizeType min_freq, SizeType max_freq, SizeType ppi) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        SizeType n = f.result_size();
        SizeType m = g.size();

        SizeType pps = n+ppi*m;

        for (auto freq : range(min_freq,max_freq+1)) {

            SizeType base = avg-n-pps*(freq-1u)/2u;

            auto sweeper = make_threshold_sweeper(1e-8);
            List<SharedPointer<InclusionIntegratorApproximation>> approximations;
            approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << freq << " (b: " << base << "): " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }

    Void run_battery_fixed_relativebase(String name,
                                ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                                RealBox real_starting_set, Real evolution_time, double step, SizeType ppi, double ratio, SizeType min_freq, SizeType max_freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        SizeType n = f.result_size();
        SizeType m = g.size();

        SizeType pps = n+ppi*m;

        for (auto freq : range(min_freq,max_freq+1)) {

            SizeType base(round(ratio*pps*freq));

            auto sweeper = make_threshold_sweeper(1e-8);
            List<SharedPointer<InclusionIntegratorApproximation>> approximations;
            approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << freq << " (b: " << base << "): " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }

    Void run_battery_each_approximation(String name,
                                        ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g,
                                        RealVector noise_levels,
                                        RealBox real_starting_set, Real evolution_time, double step, SizeType freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        auto sweeper = make_threshold_sweeper(1e-8);

        for (auto approx: range(0,5)) {

            List<SharedPointer<InclusionIntegratorApproximation>> approximations;

            SizeType ppi = 0;
            switch (approx) {
                case 0:
                    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
                    ppi = 0;
                    break;
                case 1:
                    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
                    ppi = 1;
                    break;
                case 2:
                    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
                    ppi = 2;
                    break;
                case 3:
                    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
                    ppi = 2;
                    break;
                case 4:
                    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
                    ppi = 2;
                    break;
                default:
                    break;
            }

            auto n = f.result_size();
            auto m = noise.size();
            SizeType base = n + freq/2 * n + (2*freq-1)*m - (freq-1)*ppi*m/2;
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            std::cout << approximations.at(0)->getKind() << " (reset to: " << base << ")" << std::endl;

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }



    Void run_battery_fixed_absolutebase(String name,
                                        ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g,
                                        RealVector noise_levels,
                                        RealBox real_starting_set, Real evolution_time, double step, SizeType base,
                                        SizeType min_freq, SizeType max_freq) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        for (auto freq : range(min_freq,max_freq+1)) {

            auto sweeper = make_threshold_sweeper(1e-8);
            List<SharedPointer<InclusionIntegratorApproximation>> approximations;
            approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
            ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
            ValidatedConstrainedImageSet evolve_set = range(evolve_function);

            FloatDPUpperBound total_diameter(0.0);
            auto ebb = evolve_set.bounding_box();
            for (auto i : range(ebb.size())) {
                total_diameter += ebb[i].width();
            }
            std::cout << freq << ": " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }


    Void run_battery_noreset(String name,
                                        ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g, RealVector noise_levels,
                                        RealBox real_starting_set, Real evolution_time, SizeType min_step, SizeType max_step) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;


        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        for (auto s : range(min_step,max_step+1)) {

            auto step = 1.0/(2<<s);

            std::cout << "step: " << step << std::endl;

            BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
            auto sweeper = make_threshold_sweeper(1e-8);
            List<SharedPointer<InclusionIntegratorApproximation>> approximations;
            approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
            auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=1000000, number_of_variables_to_keep=1000000);

            tms start_time, end_time;
            times(&start_time);

            List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

            times(&end_time);
            clock_t ticks = end_time.tms_utime - start_time.tms_utime;
            clock_t const hz = sysconf(_SC_CLK_TCK);

            std::cout << "    " << ticks / hz << "." << ticks % hz << "s" << std::endl;
        }
    }

    Void run_test(String name, InclusionIntegratorInterface& integrator,
                  ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, RealVector noise_levels,
                  RealBox real_starting_set, Real evolution_time) const
    {
        typedef typename ValidatedVectorFunctionModelType::NumericType NumericType; typedef typename NumericType::PrecisionType PrecisionType;
        PrecisionType prec;

        BoxDomainType noise=cast_exact_box(UpperIntervalType(-1,+1)*noise_levels);
        BoxDomainType starting_set=cast_exact_box(over_approximation(real_starting_set));

        std::cout << "flowing..." << std::endl;

        tms start_time, end_time;
        times(&start_time);

        List<ValidatedVectorFunctionModelType> flow_functions = integrator.flow(f,g,noise,starting_set,evolution_time);

        times(&end_time);
        clock_t ticks = end_time.tms_utime - start_time.tms_utime;
        clock_t const hz = sysconf(_SC_CLK_TCK);

        List<ValidatedConstrainedImageSet> reach_sets = map([](ValidatedVectorFunctionModelType const& fm){return range(fm);},flow_functions);
        ValidatedVectorFunctionModelType evolve_function = partial_evaluate(flow_functions.back(),starting_set.size(),NumericType(evolution_time,prec));
        ValidatedConstrainedImageSet evolve_set = range(evolve_function);

        FloatDPUpperBound total_diameter(0.0);
        auto ebb = evolve_set.bounding_box();
        for (auto i : range(ebb.size())) {
            total_diameter += ebb[i].width();
        }
        std::cout << "total diameter: " << total_diameter << ", " << ticks / hz << "." << ticks % hz << "s" << std::endl;

        std::cout << "plotting..." << std::endl;
        Box<FloatDPUpperInterval> graphics_box(f.result_size());
        for (auto set: reach_sets) {
            graphics_box = hull(graphics_box,set.bounding_box());
        }
        for (SizeType i : range(0,f.result_size()-1)) {
            for (SizeType j : range(i+1,f.result_size())) {
                Figure fig=Figure();
                fig.set_bounding_box(graphics_box);
                fig.set_projection(f.result_size(),i,j);
                fig.set_line_colour(0.0,0.0,0.0);
                fig.set_line_style(false);
                fig.set_fill_colour(0.5,0.5,0.5);
                fig.draw(starting_set);
                fig.set_fill_colour(1.0,0.75,0.5);
                for (auto set : reverse(reach_sets)) { fig.draw(set); }
                fig.draw(evolve_set);
                char num_char[7] = "";
                if (f.result_size() > 2)
                    sprintf(num_char,"[%lu,%lu]",i,j);
                fig.write(("test_differential_inclusion-"+name+num_char).c_str());
            }
        }
    }

  public:
    void test() const;

    void test_reactor() const;
    void test_lorenz() const;
    void test_rossler() const;
    void test_jerk21() const;
    void test_jerk16() const;
    void test_higgins_selkov() const;
    void test_lotka_volterra() const;
    void test_jet_engine() const;
    void test_pi_controller() const;
    void test_van_der_pol() const;
    void test_harmonic() const;
    void test_harmonic_analytical() const;
    void test_clock() const;
    void test_DCDC() const;
};

void TestInclusionIntegrator::test_reactor() const {
    double step=1.0/32;
    SizeType freq=12;
    SizeType base=1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/100_q,1/100_q,1/100_q};

    auto x = EffectiveVectorFunction::identity(4u);
    auto one = EffectiveScalarFunction::constant(4u,1_z);
    auto zero = EffectiveScalarFunction::constant(4u,0_z);

    Real U3(30.0);
    Real k2(0.4);
    Real iV(0.05);
    Real ka(0.1);
    Real kb(0.045);

    auto f = EffectiveVectorFunction({-U3*x[0]*x[1]-k2*x[0]*x[2]+iV-ka*x[0],-U3*x[0]*x[1]+kb-ka*x[1],
                                      U3*x[0]*x[1]-k2*x[0]*x[2]-ka*x[2],k2*x[0]*x[2]-ka*x[3]});

    Vector<ValidatedVectorFunction> g({{one*iV,zero,zero,zero},{zero,one*iV,zero,zero},{x[0]*x[1],x[0]*x[1],x[0]*x[1],zero}});

    Real e=1/1000000_q;
    RealBox starting_set={{0,e},{0,e},{0,e},{0,e}};
    Real evolution_time=80/10_q;

    this->run_test("reactor",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("reactor",f,g,noise_levels,starting_set,evolution_time,step,freq);
}


void TestInclusionIntegrator::test_jerk16() const {
    double step=1.0/16;
    SizeType freq=12;
    SizeType n=3;
    SizeType m=1;
    SizeType params_per_step=2;
    SizeType base = n + freq/2 * n + (2*freq-1)*m - (freq-1)*params_per_step*m/2;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1000_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real B(0.03);

    auto f = EffectiveVectorFunction({x[1],x[2],-x[1]+x[0]*x[0]-one*B});

    Vector<ValidatedVectorFunction> g({{zero,zero,one}});

    Real e=1/1024_q;
    Real x0_i(0.0);
    Real x1_i(0.0);
    Real x2_i(0.0);
    RealBox starting_set={{x0_i+e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=100/10_q;

    this->run_test("jerk16",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("jerk16",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_jerk21() const {

    double step=1.0/16;
    SizeType freq=12;
    SizeType base=1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1000_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real A(0.25);

    auto f = EffectiveVectorFunction({x[1],x[2],-x[2]*x[2]*x[2]-x[1]*x[0]*x[0]-A*x[0]});

    Vector<ValidatedVectorFunction> g({{zero,zero,-x[0]}});

    Real e=1/1024_q;
    Real x0_i(0.25);
    Real x1_i(0.0);
    Real x2_i(0.0);
    RealBox starting_set={{x0_i+e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=100/10_q;

    this->run_test("jerk21",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("jerk21",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_lorenz() const {
    double step=1.0/256;
    SizeType freq=12;
    SizeType base = 1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/100_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real sigma(10.0);
    Real rho(28.0);
    Real beta(8.0/3);

    auto f = EffectiveVectorFunction({(x[1]-x[0])*sigma,x[0]*(one*rho - x[2]) - x[1],x[0]*x[1] - x[2]*beta});

    Vector<ValidatedVectorFunction> g({{zero,x[0],zero}});

    Real e=1/1024_q;
    Real x0_i(1.0);
    Real x1_i(1.0);
    Real x2_i(1.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=10/10_q;

    this->run_test("lorenz",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("lorenz",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_rossler() const {
    double step=1.0/128;
    SizeType freq=12;
    SizeType n=3;
    SizeType m=1;
    SizeType params_per_step=2;
    SizeType base = n + freq/2 * n + (2*freq-1)*m - (freq-1)*params_per_step*m/2;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/1000_q};

    auto x = EffectiveVectorFunction::identity(3u);
    auto one = EffectiveScalarFunction::constant(3u,1_z);
    auto zero = EffectiveScalarFunction::constant(3u,0_z);

    Real a(0.1);
    Real b(0.1);
    Real c(6.0);

    auto f = EffectiveVectorFunction({-x[1]-x[2],x[0] + x[1]*a,one*b + x[2]*(x[0]-one*c)});

    Vector<ValidatedVectorFunction> g({{zero,zero,one}});

    Real e=1/1024_q;
    Real x0_i(-9.0);
    Real x1_i(0.0);
    Real x2_i(0.01);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e},{x2_i-e,x2_i+e}};
    Real evolution_time=120/10_q;

    this->run_test("rossler",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("rossler",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_higgins_selkov() const {

    double step=1.0/50;
    SizeType freq=12;
    SizeType base=1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={2/10000_q,2/10000_q,2/10000_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    Real k1(1.00001);

    auto f = EffectiveVectorFunction({one-x[0]*k1*x[1]*x[1],x[0]*k1*x[1]*x[1] - x[1]});

    Vector<ValidatedVectorFunction> g({{one,zero},{-x[0]*x[1]*x[1],x[0]*x[1]*x[1]},{zero,-x[1]}});

    Real e=1/100_q;
    Real x0_i(2.0);
    Real x1_i(1.0);
    RealBox starting_set={{x0_i-e,x0_i+e},{x1_i-e,x1_i+e}};
    Real evolution_time=100/10_q;

    this->run_test("higgins-selkov",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("higgins-selkov",f,g,noise_levels,starting_set,evolution_time,step,freq);
}


void TestInclusionIntegrator::test_lotka_volterra() const {
    double step=1.0/50;
    SizeType freq=12;
    SizeType base=1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/100_q,1/100_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    auto three = EffectiveScalarFunction::constant(2u,3_z);

    auto f = EffectiveVectorFunction({three*x[0]*(one-x[1]),x[1]*(x[0]-one)});

    Vector<ValidatedVectorFunction> g({{x[0]*(one-x[1]),zero},{zero,x[1]*(x[0]-one)}});

    Real e=1/100000000_q;
    RealBox starting_set={{Real(1.2)-e,Real(1.2)+e},{Real(1.1)-e,Real(1.1)+e}};
    Real evolution_time=100/10_q;

    this->run_test("lotka-volterra",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("lotka-volterra",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_jet_engine() const {
    double step=1.0/50;
    SizeType freq=12;
    SizeType base = 1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={5/1000_q,5/1000_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f = EffectiveVectorFunction({-x[1]-Real(1.5)*x[0]*x[0]-Real(0.5)*x[0]*x[0]*x[0]-Real(0.5),3*x[0]-x[1]});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e1=5/100_q;
    Real e2=7/100_q;
    RealBox starting_set={{Real(1.0)-e1,Real(1.0)+e1},{Real(1.0)-e2,Real(1.0)+e2}};
    Real evolution_time=40/8_q;

    this->run_test("jet-engine",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("jet-engine",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_pi_controller() const {
    double step=1.0/32;
    SizeType freq=12;
    SizeType base = 1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/10_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f = EffectiveVectorFunction({Real(-0.101)*(x[0]-Real(20.0))+Real(1.3203)*(x[1]-Real(0.1616))-Real(0.01)*x[0]*x[0], Real(-1.0)*(Real(-0.101)*(x[0]-Real(20.0))+Real(1.3203)*(x[1]-Real(0.1616))-Real(0.01)*x[0]*x[0]) + Real(3.0)*(Real(20.0)-x[0])});

    Vector<ValidatedVectorFunction> g({{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(5.0),Real(10.0)},{-e,+e}};
    Real evolution_time=40/8_q;

    this->run_test("pi-controller",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("pi-controller",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_van_der_pol() const {

    double step=1.0/8;
    SizeType freq=10;
    SizeType base=16;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/20_q,1/10000_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({x[1],-x[0]+x[1]*(1 - x[0]*x[0])});

    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/1024_q;
    RealBox starting_set={{Real(1.21)-e,Real(1.21)+e},{Real(2.01)-e,Real(2.01)+e}};
    Real evolution_time=8/4_q;

    this->run_test("vanderpol",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test_DCDC() const {
    double step=1.0/10;
    SizeType freq=12;
    SizeType base = 1000;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorPiecewiseApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorSinusoidalApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorConstantApproximation(sweeper)));
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorZeroApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={2/1000_q,1/15_q};

    Real k0(0.002987);
    Real fp0(-0.018);//(-11+k0)/600
    Real fp1(-0.066); //(k0-1)/15
    Real fq0(0.071);//(1-k0)/14
    Real fq1(-0.00853);//-k0*20/7
    Real gp0 = 1/600_q;
    Real gp1 = 1/15_q;
    Real gq0 = -1/14_q;
    Real gq1 = -20/7_q;

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);
    auto third = EffectiveScalarFunction::constant(2u,1/3_z);

    auto f = EffectiveVectorFunction({third+x[0]*fp0+x[1]*fp1,x[0]*fq0+x[1]*fq1});

    Vector<ValidatedVectorFunction> g({{gp0*x[0]+gp1*x[1],gq0*x[0]+gq1*x[1]},{one,zero}});

    Real e=1/1000000_q;
    RealBox starting_set={{Real(1)-e,Real(1)+e},{Real(5)-e,Real(5)+e}};
    Real evolution_time=50/10_q;

    this->run_test("DCDC",integrator,f,g,noise_levels,starting_set,evolution_time);
    //this->run_battery_each_approximation("DCDC",f,g,noise_levels,starting_set,evolution_time,step,freq);
}

void TestInclusionIntegrator::test_harmonic() const {
    double step=1.0/64;
    SizeType base=80;
    SizeType freq=1000;

    SweeperDP sweeper=make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 0;

    RealVector noise_levels={4/100_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({x[1],-x[0]});

    Vector<ValidatedVectorFunction> g({{one,zero}});

    Real e=1/10000000_q;
    RealBox starting_set={{-e,e},{-e,+e}};
    Real evolution_time=314/100_q;

    //this->run_test("harmonic",integrator,f,g,noise_levels,starting_set,evolution_time);
    this->run_battery_each_approximation("harmonic",f,g,noise_levels,starting_set,evolution_time,step,freq);
}


void TestInclusionIntegrator::test_clock() const {

    double step=1.0/256;
    SizeType base=24;
    SizeType freq=11;

    auto sweeper = make_threshold_sweeper(1e-8);
    List<SharedPointer<InclusionIntegratorApproximation>> approximations;
    approximations.append(SharedPointer<InclusionIntegratorApproximation>(new InclusionIntegratorAffineApproximation(sweeper)));
    auto integrator = InclusionIntegrator(approximations,sweeper,step_size=step, number_of_steps_between_simplifications=freq, number_of_variables_to_keep=base);
    integrator.verbosity = 2;

    RealVector noise_levels={1/16_q,1/16_q};

    auto x = EffectiveVectorFunction::identity(2u);
    auto one = EffectiveScalarFunction::constant(2u,1_z);
    auto zero = EffectiveScalarFunction::constant(2u,0_z);

    auto f=EffectiveVectorFunction({one,one});
    Vector<ValidatedVectorFunction> g({{one,zero},{zero,one}});

    Real e=1/128_q;
    RealBox starting_set={{-e,e},{-e,+e}};
    Real evolution_time=20/4_q;

    this->run_test("clock",integrator,f,g,noise_levels,starting_set,evolution_time);
}

void TestInclusionIntegrator::test() const {
    //ARIADNE_TEST_CALL(test_reactor());
    //ARIADNE_TEST_CALL(test_jerk16());
    //ARIADNE_TEST_CALL(test_jerk21());
    //ARIADNE_TEST_CALL(test_lorenz());
    //ARIADNE_TEST_CALL(test_rossler());
    ARIADNE_TEST_CALL(test_higgins_selkov());
    //ARIADNE_TEST_CALL(test_lotka_volterra());
    //ARIADNE_TEST_CALL(test_jet_engine());
    //ARIADNE_TEST_CALL(test_pi_controller());
    //ARIADNE_TEST_CALL(test_van_der_pol());
    //ARIADNE_TEST_CALL(test_DCDC());
    //ARIADNE_TEST_CALL(test_harmonic());
    //ARIADNE_TEST_CALL(test_clock());
}

int main() {
    TestInclusionIntegrator().test();
    return ARIADNE_TEST_FAILURES;
}
