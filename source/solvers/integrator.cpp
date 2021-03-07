/***************************************************************************
 *            solvers/integrator.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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

#include "function/functional.hpp"
#include "config.hpp"

#include <iomanip>

#include "solvers/integrator.hpp"
#include "solvers/bounder.hpp"

#include "output/logging.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/differential.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"
#include "function/formula.hpp"
#include "function/scaling.hpp"
#include "function/taylor_function.hpp"
#include "function/taylor_model.hpp"

#include "function/polynomial.hpp"
#include "geometry/interval.hpp"

#include "algebra/expansion.inl.hpp"

namespace Ariadne {

typedef ValidatedVectorMultivariateTaylorFunctionModelDP FlowStepTaylorModelType;


OutputStream& operator<<(OutputStream& os, FlowStepModelType const& fsm) {
    return os << static_cast<ValidatedVectorMultivariateFunctionModelDP const&>(fsm);
}

OutputStream& operator<<(OutputStream& os, FlowModelType const& fm) {
    os << "[ "; for(SizeType i=0; i!=fm.size(); ++i) { if (i!=0u) { os << ",\n"; } os << "  " << fm[i]; } os << " ]"; return os ;
}

static const FloatDPValue zero=FloatDPValue(0,dp);

inline UpperBoxType operator+(Vector<ExactIntervalType> bx, Vector<UpperIntervalType> const& ex) {
    return Vector<UpperIntervalType>(bx) + ex;
}

inline UpperBoxType operator+(Vector<UpperIntervalType> bx, Vector<FloatDPBounds> const& v) {
    return bx + Vector<UpperIntervalType>(v);
}

inline UpperBoxType operator+(Vector<ExactIntervalType> bx, Vector<FloatDPBounds> const& v) {
    return Vector<UpperIntervalType>(bx) + Vector<UpperIntervalType>(v);
}

inline ExactDouble cast_exact_double(Attribute<ApproximateDouble> a) { return cast_exact(static_cast<ApproximateDouble>(a)); }

IntegratorBase::IntegratorBase(Configuration<IntegratorBase> const& config)
    : Configurable<IntegratorBase>(config) { }

StepSizeType
IntegratorBase::starting_time_step_size(ValidatedVectorMultivariateFunction const& vector_field, BoxDomainType const& state_domain, StepSizeType const& evaluation_time_step, StepSizeType const& maximum_time_step) const {
    FloatDPUpperBound lipschitz = norm(vector_field.jacobian(Vector<FloatDPBounds>(cast_singleton(product(state_domain, to_time_bounds(0, evaluation_time_step)))))).upper();
    auto lipschitz_step = static_cast<StepSizeType>(cast_exact(configuration().lipschitz_tolerance()/lipschitz));
    ARIADNE_LOG_PRINTLN_VAR_AT(1,lipschitz_step);
    return min(lipschitz_step, maximum_time_step);
}

StepSizeType
IntegratorBase::starting_time_step_size(ValidatedVectorMultivariateFunction const& vector_field, BoxDomainType const& state_domain, BoxDomainType const& parameter_domain, StepSizeType const& evaluation_time_step, StepSizeType const& maximum_time_step) const {

    FloatDPUpperBound lipschitz = norm(vector_field.jacobian(Vector<FloatDPBounds>(cast_singleton(product(state_domain, to_time_bounds(0, evaluation_time_step), parameter_domain))))).upper();
    return min(static_cast<StepSizeType>(cast_exact(configuration().lipschitz_tolerance()/lipschitz)), maximum_time_step);
}

Pair<StepSizeType,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& hsug) const {
    return configuration().bounder().compute(vf,D,hsug);
}

Pair<StepSizeType,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const ExactBoxType& A, const StepSizeType& hsug) const {
    return configuration().bounder().compute(vf,D,A,hsug);
}

Pair<StepSizeType,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, StepSizeType const& t, const ExactBoxType& A, const StepSizeType& hsug) const {
    return configuration().bounder().compute(vf,D,t,A,hsug);
}


FlowStepModelType
IntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& hev, StepSizeType& hsug) const
{
    StepSizeType& h=hsug;
    h = starting_time_step_size(vf, dx, hev, h);
    UpperBoxType bx;
    make_lpair(h,bx)=this->flow_bounds(vf,dx,h);
    while(true) {
        try {
            return this->flow_step(vf,dx,h,bx);
        } catch(const FlowTimeStepException& e) {
            h=hlf(h);
            ARIADNE_LOG_PRINTLN("Failure to compute the flow step for the current step size, halving it to " << h);
        }
    }
}

inline ExactDouble operator*(ExactDouble, TwoExp);
inline ExactDouble operator/(ExactDouble, TwoExp);

TaylorPicardIntegrator::TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator> const& config)
    : IntegratorBase(config) { }

TaylorPicardIntegrator* TaylorPicardIntegrator::clone() const {
    return new TaylorPicardIntegrator(this->configuration());
}

Configuration<TaylorPicardIntegrator> const& TaylorPicardIntegrator::configuration() const {
    return static_cast<Configuration<TaylorPicardIntegrator> const&>(IntegratorBase::configuration());
}

FlowStepModelType
TaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& h, const UpperBoxType& bx) const
{
    ARIADNE_PRECONDITION(vf.result_size()==dx.dimension());
    ARIADNE_PRECONDITION(vf.argument_size()==dx.dimension());
    ARIADNE_PRECONDITION(bx.dimension()==dx.dimension());
    return this->_flow_step(vf,dx,IntervalDomainType(0,h),BoxDomainType(0u),bx);
}

FlowStepModelType
TaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+T.dimension()+A.dimension());
    ARIADNE_PRECONDITION(B.dimension()==D.dimension());
    return this->_flow_step(f,D,IntervalDomainType(T),A,B);
}

FlowStepModelType
TaylorPicardIntegrator::_flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const ExactIntervalType& T, const ExactBoxType& A, const UpperBoxType& B) const {
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("f=" << f);
    ARIADNE_LOG_PRINTLN("D=" << D << " T=" << T << ", A=" << A << ", B=" << B);

    const bool is_autonomous = (f.argument_size() == D.dimension() + A.dimension());

    const SizeType nx = D.size();
    const SizeType na = A.size();

    Range tarng = is_autonomous ? Range(nx + 1u, nx + 1u + na) : Range(nx, nx + 1u + na);

    StepSizeType t = static_cast<StepSizeType>(T.lower_bound());
    StepSizeType h = static_cast<StepSizeType>(T.upper_bound()) - t;

    // Time interval centred on initial time, which will make the antiderivative more efficient
    ExactIntervalType wT(t - h, t + h);
    ARIADNE_ASSERT(t == med(wT));

    ExactBoxType dom = join(D, T, A);
    ExactBoxType wdom = join(D, wT, A);
    UpperBoxType const &bx = B;
    ARIADNE_LOG_PRINTLN_AT(2, "dom=" << dom << ", wdom=" << wdom);

    auto function_factory = make_taylor_function_factory(configuration().sweeper());

    FlowStepModelType phi0 = function_factory->create_projection(wdom, range(0, nx));
    ARIADNE_LOG_PRINTLN_VAR_AT(1, phi0);
    FlowStepModelType phi = function_factory->create_constants(wdom, cast_singleton(bx));
    FlowStepModelType ta = function_factory->create_projection(wdom, tarng);

    ARIADNE_LOG_PRINTLN_VAR_AT(1, phi);
    for (DegreeType k = 0; k != configuration().maximum_temporal_order(); ++k) {
        Bool below_maximum_error = (phi.error().raw() < configuration().step_maximum_error());

        FlowStepModelType fphi = compose(f, join(std::move(phi), ta));
        ARIADNE_LOG_PRINTLN_VAR_AT(2, fphi);
        // NOTE: In principle safer to use antiderivative(fphi,nx,t) here,
        // but since t is the midpoint of wdom, the (standard) antiderivative works
        // TODO: Change based antiderivative to be efficient when t is midpoint of domain
        phi = antiderivative(fphi, nx) + phi0;
        ARIADNE_LOG_PRINTLN_VAR_AT(2, phi);
        if (below_maximum_error && k >= configuration().minimum_temporal_order()) { break; }
    }
    if (phi.error().raw() > configuration().step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException, "TaylorPicardIntegrator::flow_step",
                      "Integration of " << f << " starting in " << D << " over time interval " << T << " of length "
                                        << h
                                        << " has error " << phi.error() << " after "
                                        << configuration().maximum_temporal_order() <<
                                        " iterations, which exceeds step maximum error "
                                        << configuration().step_maximum_error());
    }

    FlowStepModelType result = restriction(phi, dom);

    ARIADNE_LOG_PRINTLN_VAR(result);
    return result;
}


Void TaylorPicardIntegrator::_write(OutputStream& os) const {
    os << "TaylorPicardIntegrator" << configuration();
}

} // namespace Ariadne


#include "algebra/graded.hpp"
#include "function/procedure.hpp"
namespace Ariadne {

class FormulaFunction;
typedef Procedure<ValidatedNumber> ValidatedProcedure;
typedef Differential<FloatDPBounds> ValidatedDifferential;
typedef MultivariatePolynomial<FloatDPBounds> ValidatedMultivariatelynomial;
typedef Graded<ValidatedDifferential> GradedValidatedDifferential;
typedef FloatDPBounds ValidatedNumericType;
Bool operator<(const MultiIndex& a1, const MultiIndex& a2);

static const TwoExp STEP_ERROR_RATIO=TwoExp(-7);
static const TwoExp SWEEP_THRESHOLD_RATIO=TwoExp(-10);
static const LipschitzConstant DEFAULT_LIPSCHITZ_CONSTANT(0.5_x);

TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip, Order ord)
    : IntegratorBase(/* FIXME */ Configuration<TaylorPicardIntegrator>()), _sweeper(sweeper), _order(ord)
{ }

TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        MaximumError err, Order ord)
    : TaylorSeriesIntegrator(err,ThresholdSweeper<FloatDP>(DP(),err.value()*SWEEP_THRESHOLD_RATIO),DEFAULT_LIPSCHITZ_CONSTANT,ord)
{ }


GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(MaximumError err)
    : GradedTaylorSeriesIntegrator(err,ThresholdSweeper<FloatDP>(DP(),err.value()*SWEEP_THRESHOLD_RATIO),DEFAULT_LIPSCHITZ_CONSTANT)
{ }

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip)
    : GradedTaylorSeriesIntegrator(err,sweeper,lip,StepMaximumError(err.value()*STEP_ERROR_RATIO),MaximumTemporalOrder(12))
{ }

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                        StepMaximumError stperr, MaximumTemporalOrder maxto)
    : GradedTaylorSeriesIntegrator(err,sweeper,lip,stperr,MinimumSpacialOrder(1),MinimumTemporalOrder(4),MaximumSpacialOrder(4),maxto)
{
}

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(
        MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip, StepMaximumError stperr,
        MinimumSpacialOrder minso, MinimumTemporalOrder minto,
        MaximumSpacialOrder maxso, MaximumTemporalOrder maxto)
    : IntegratorBase(/* FIXME */Configuration<TaylorPicardIntegrator>()), _step_maximum_error(cast_exact(stperr.value())), _sweeper(sweeper)
    , _minimum_spacial_order(minso), _minimum_temporal_order(minto), _maximum_spacial_order(maxso), _maximum_temporal_order(maxto)
{ }


namespace {

ExactIntervalType forwards_backwards_time_domain(ExactIntervalType domt) {
    Dyadic t0(domt.lower_bound());
    Dyadic tf(domt.upper_bound());
    return ExactIntervalType(t0-(tf-t0),tf);
}

template<class F> GradedValidatedDifferential flow(const F& f, const ExactIntervalType& c, DegreeType so, DegreeType to) {
    ValidatedDifferential x=make_differential_variable(1u,so,cast_singleton(c),0u);
    GradedValidatedDifferential y=make_graded(x);
    GradedValidatedDifferential t=create_graded(x);

    for(DegreeType n=0; n!=to; ++n) {
        t=f(y);
        y=antidifferential(t);
    }

    return y;
}

template<class X> Void append_join(Expansion<MultiIndex,X>& e, const MultiIndex& a1, const DegreeType a2, const X& c) {
    MultiIndex a(a1.size()+1);
    for(SizeType i=0; i!=a1.size(); ++i) { a[i]=a1[i]; }
    a[a1.size()]=a2;
    e.append(a,c);
}



Void graded_flow_init(const Vector<ValidatedProcedure>& f,
               Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp,
               Vector<GradedValidatedDifferential>& yta,
               const Vector<ValidatedNumericType>& x, const ValidatedNumericType& t0, const Vector<ValidatedNumericType>& a,
               DegreeType so, DegreeType to)
{
    ARIADNE_LOG_SCOPE_CREATE;
    const SizeType xs=x.size();
    const SizeType as=a.size();
    const SizeType ress=f.result_size();
    const SizeType tmps=f.temporaries_size();
    const SizeType args=f.argument_size();

    assert(args==xs+1u+as or args==xs+as);
    const bool is_autonomous = (args==xs+as);

    GradedValidatedDifferential null;
    fy=Vector< GradedValidatedDifferential >(ress,null);
    tmp=List< GradedValidatedDifferential >(tmps,null);
    yta=Vector< GradedValidatedDifferential >(args,null);
    for(SizeType i=0; i!=xs; ++i) {
        yta[i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,x[i],i));
    }
    if (is_autonomous) {
        for(SizeType i=0; i!=as; ++i) {
            yta[xs+i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,a[i],xs+i));
        }
    } else {
        yta[xs]=GradedValidatedDifferential(ValidatedDifferential::constant(xs+as,so,t0));
        for(SizeType i=0; i!=as; ++i) {
            yta[xs+1u+i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,a[i],xs+i));
        }
    }
    ARIADNE_LOG_PRINTLN_AT(1,"fy="<<fy<<", tmp="<<tmp<<", yta="<<yta);
}


Void graded_flow_iterate(const Vector<ValidatedProcedure>& p,
                         Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp, Vector<GradedValidatedDifferential>& yta)
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"degree="<<yta[0].degree());
    const bool is_autonomous = (p.argument_size()==yta[0][0].argument_size());
    const SizeType n=p.result_size();

    ValidatedDifferential z=nul(yta[0][0]);

    Ariadne::compute(p,fy,tmp,yta);
    for(SizeType i=0; i!=n; ++i) {
        yta[i]=antidifferential(fy[i]);
    }

    if (is_autonomous) {
        for(SizeType i=n; i!=yta.size(); ++i) { yta[i].append(z); }
    } else {
        GradedValidatedDifferential& t=yta[n];
        if(t.degree()==0 && not is_autonomous) { t.append(z+1); } else { t.append(z); }
        for(SizeType i=n+1u; i!=yta.size(); ++i) { yta[i].append(z); }
    }
}




Vector<GradedValidatedDifferential>
graded_flow_differential(Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphib,
                         DegreeType so, DegreeType to)
{
    const SizeType rs=dphic.size();
    const SizeType as=dphic[0][0].argument_size();
    auto z=dphic[0][0].zero_coefficient();

    Vector<GradedValidatedDifferential> gdphi(rs,GradedValidatedDifferential(List<ValidatedDifferential>(to+1u,ValidatedDifferential(as,so,z))));
    for(SizeType i=0; i!=rs; ++i) {
        for(DegreeType j=0; j!=to; ++j) {
            for(ValidatedDifferential::ConstIterator iter=dphic[i][j].begin(); iter!=dphic[i][j].end(); ++iter) {
                if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
            for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
                if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
        }
        DegreeType j=to;
        for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
            gdphi[i][j].expansion().append(iter->index(),iter->coefficient());
        }
    }
    ARIADNE_LOG_PRINTLN_AT(1,"gdphi="<<gdphi);

    return gdphi;
}

Vector<ValidatedDifferential>
differential(Vector<GradedValidatedDifferential> const& gdphi, SizeType gind,
             DegreeType so, DegreeType to)
{
    SizeType rs=gdphi.size();
    SizeType as=gdphi[0][0].argument_size();
    auto z=gdphi[0][0].zero_coefficient();

    Vector<ValidatedDifferential> dphi(rs,as+1u,so+to,z);
    MultiIndex a(as+1u);
    for(SizeType i=0; i!=rs; ++i) {
        Expansion<MultiIndex,FloatDPBounds>& component=dphi[i].expansion();
        for(DegreeType j=0; j<=to; ++j) {
            a[gind]=j;
            const Expansion<MultiIndex,FloatDPBounds>& expansion=gdphi[i][j].expansion();
            for(auto term : expansion) {
                for(SizeType k=0; k!=gind; ++k) { a[k]=term.index()[k]; }
                for(SizeType k=gind; k!=as; ++k) { a[k+1u]=term.index()[k]; }
                component.append(a,term.coefficient());
            }
        }
    }
    ARIADNE_LOG_PRINTLN_AT(1,"dphi="<<dphi);
    return dphi;
}

Vector<ValidatedDifferential>
flow_differential(Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphib,
                  DegreeType so, DegreeType to)
{
    Vector<GradedValidatedDifferential> gdphi=graded_flow_differential(dphic,dphib,so,to);
    return differential(gdphi, dphic.size(),so,to);
}


FlowStepTaylorModelType make_taylor_function_model(const Vector<Differential<FloatBounds<DP>>>& df, const ExactBoxType& dom, Sweeper<FloatDP> swp) {
    ARIADNE_ASSERT(df.argument_size()==dom.dimension());
    const SizeType rs=df.size();
    const SizeType as=dom.dimension();
    const DegreeType deg = df.degree();
    FlowStepTaylorModelType tf(rs,dom,swp);

    Vector<Differential<FloatBounds<DP>>> ds=scale(Differential<FloatBounds<DP>>::variables(deg,Vector<FloatBounds<DP>>(as,dp)),dom);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,ds);
    Vector<Differential<FloatBounds<DP>>> dfs = compose(df,ds);

    for(SizeType i=0; i!=rs; ++i) {
        ValidatedTaylorModelDP& model=tf.model(i);
        Expansion<MultiIndex,FloatDPValue>& expansion=model.expansion();
        FloatDPError& error=model.error();
        error=0u;
        expansion.reserve(dfs[i].expansion().number_of_nonzeros());

        typename Differential<FloatDPBounds>::ConstIterator iter=dfs[i].begin();
        while(iter!=dfs[i].end()) {
            MultiIndex const a=iter->index();
            FloatDPBounds coef=iter->coefficient();
            FloatDPValue x=coef.value();
            error+=coef.error();
            expansion.append(a,x);
            ++iter;
        }
        model.cleanup();
    }
    return tf;
}

FlowStepTaylorModelType flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, Sweeper<FloatDP> swp) {
    StepSizeType t=static_cast<StepSizeType>(domt.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(domt.upper_bound())-t;
    ExactIntervalType wdt(t-h,t+h);

    return restriction(make_taylor_function_model(dphi,join(domx,wdt,doma),swp),join(domx,domt,doma));
}

} // namespace


// Flow step using graded differential with fixed degree
FlowStepTaylorModelType
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        Sweeper<FloatDP> const& sweeper, DegreeType so, DegreeType to)
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"f="<<f);
    ARIADNE_LOG_PRINTLN_AT(1,"domx="<<domx<<", domt="<<domt<<", doma="<<doma<<", bndx="<<bndx);
    ARIADNE_LOG_PRINTLN_AT(1,"sweeper="<<sweeper<<", so="<<so<<", to="<<to);

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension() || f.argument_size()==domx.dimension()+1u+doma.dimension());

    SizeType nx=domx.dimension();

    StepSizeType t=static_cast<StepSizeType>(domt.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(domt.upper_bound())-t;

    ExactIntervalType widt(t-h,t+h);

    Vector<ValidatedNumericType> dx=cast_singleton(domx);
    Scalar<ValidatedNumericType> dt=cast_singleton(domt);
    Scalar<ValidatedNumericType> wdt=cast_singleton(widt);
    Vector<ValidatedNumericType> da=cast_singleton(doma);
    Vector<ValidatedNumericType> bx=cast_singleton(bndx);

    Vector<ValidatedNumericType> mdx=midpoint(dx);
    Scalar<ValidatedNumericType> mdt=midpoint(widt);
    Vector<ValidatedNumericType> mda=midpoint(da);

    Vector<ValidatedNumericType> dc=midpoint(join(dx,wdt,da));
    Vector<ValidatedNumericType> db=join(bx,dt,da);
    ExactBoxType domc=join(domx,domt,doma);

    ARIADNE_LOG_PRINTLN_AT(2,"dx="<<dx<<", dt="<<dt<<", da="<<da<<", wdt="<<wdt<<", bx="<<bx);

    Vector<GradedValidatedDifferential> dphic,fdphic,dphib,fdphib;
    List<GradedValidatedDifferential> tmpdphic,tmpdphib;

    Ariadne::graded_flow_init(f,fdphic,tmpdphic,dphic,mdx,mdt,mda,so,to);
    Ariadne::graded_flow_init(f,fdphib,tmpdphib,dphib,bx,dt,da,so,to);

    for(DegreeType i=0; i!=to; ++i) {
        Ariadne::graded_flow_iterate(f,fdphic,tmpdphic,dphic);
        Ariadne::graded_flow_iterate(f,fdphib,tmpdphib,dphib);
    }
    ARIADNE_LOG_PRINTLN_AT(3,"dphic="<<dphic);
    ARIADNE_LOG_PRINTLN_AT(3,"dphib="<<dphib);

    dphic=project(dphic,range(0,nx));
    dphib=project(dphib,range(0,nx));
    ARIADNE_LOG_PRINTLN_AT(3,"dphic="<<dphic);
    ARIADNE_LOG_PRINTLN_AT(3,"dphib="<<dphib);

    Vector<ValidatedDifferential> dphi=Ariadne::flow_differential(dphic,dphib,so,to);
    ARIADNE_LOG_PRINTLN_AT(2,"dphi="<<dphi);

    FlowStepTaylorModelType tphi=Ariadne::flow_function(dphi,domx,domt,doma,sweeper);

    ARIADNE_LOG_PRINTLN("phi="<<tphi);

    return tphi;
}

// Flow step using graded differential with varying degree and specified maximum error
FlowStepModelType
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        ExactDouble max_err, Sweeper<FloatDP> const& sweeper, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to)
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"f="<<f);
    ARIADNE_LOG_PRINTLN_AT(1,"domx="<<domx<<", domt="<<domt<<", doma="<<doma<<", bndx="<<bndx);
    ARIADNE_LOG_PRINTLN_AT(1,"max_err="<<max_err<<", sweeper="<<sweeper<<", "<<
                            "init_so="<<init_so<<", init_to="<<init_to<<", max_so="<<max_so<<", max_to="<<max_to);

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension() || f.argument_size()==domx.dimension()+1u+doma.dimension());

    static const ExactDouble TRY_SPACIAL_ORDER_INCREASE_FACTOR=4;

    DegreeType so=init_so;
    DegreeType to=init_to;

    FlowStepTaylorModelType phi=graded_series_flow_step(f,domx,domt,doma,bndx, sweeper,so,to);

    ARIADNE_LOG_PRINTLN_AT(1,"phi="<<phi);
    SizeType nnz=0; for(SizeType i=0; i!=phi.size(); ++i) { nnz+=phi.model(i).number_of_nonzeros(); }
    ARIADNE_LOG_PRINTLN_AT(1,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<phi.error());

    FloatDPError old_error=phi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR*two,dp);

    while(phi.error().raw()>max_err && (so<max_so || to<max_to) ) {

        old_error = phi.error();

        if( (so<max_so) && ((phi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR,dp)).raw() > old_error.raw()) ) {
            // try increasing spacial degree
            ++so;
            to=init_to;
        } else {
            ++to;
        }

        phi=graded_series_flow_step(f,domx,domt,doma,bndx, sweeper,so,to);

        SizeType nnnz=0; for(SizeType i=0; i!=phi.size(); ++i) { nnnz+=phi.model(i).number_of_nonzeros(); }
        ARIADNE_LOG_PRINTLN_AT(2,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<phi.error());
    }
    ARIADNE_LOG_PRINTLN("phi="<<phi);
    return phi;
}


FlowStepModelType
graded_series_flow_step(const ValidatedVectorMultivariateFunction& f,
                        const ExactBoxType& domx, const Interval<StepSizeType>& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        ExactDouble max_err, Sweeper<FloatDP> const& sweeper, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to)
{
    Vector<ValidatedProcedure> p(f);
    ExactIntervalType idomt(domt);

    return graded_series_flow_step(p,domx,idomt,doma,bndx,max_err,sweeper,init_so,init_to, max_so, max_to);
}


// FIXME: Should not be necessary, as should be able to construct FloatBounds<DP> from (FloatValue<DP>,DP)
FloatBounds<DoublePrecision> cast_singleton(ExactIntervalType const& ivl, DoublePrecision pr) {
    return FloatBounds<DoublePrecision>(ivl.lower_bound(),ivl.upper_bound()); }


namespace {
// Compute the midpoint of x, and add the error to e
template<class F, class FE> Value<F> med(Bounds<F> const& x, Error<FE>& e) {
    e+=x.error(); return x.value(); }
} // namespace


template<class FLT> ValidatedVectorMultivariateTaylorFunctionModel<FLT>
make_taylor_function_model(ExactBoxType domain, Vector<Differential<Bounds<FLT>>> centre_derivatives, Vector<Differential<Bounds<FLT>>> derivative_ranges, Sweeper<FLT> swp) {
    ARIADNE_ASSERT(centre_derivatives.result_size()==derivative_ranges.result_size());
    ARIADNE_ASSERT(centre_derivatives.argument_size()==domain.dimension());
    ARIADNE_ASSERT(derivative_ranges.argument_size()==domain.dimension());
    ARIADNE_ASSERT(centre_derivatives.degree()==derivative_ranges.degree() or centre_derivatives.degree()+1u==derivative_ranges.degree());

    using PR = PrecisionType<FLT>;
    using X = FloatBounds<PR>;
    PR pr = swp.precision();

    const SizeType m=derivative_ranges.result_size();
    const SizeType n=derivative_ranges.argument_size();
    const DegreeType deg = derivative_ranges.degree();
    FloatBounds<PR> z=derivative_ranges.zero_element().zero_coefficient();

    auto scalings = Vector<Differential<FloatBounds<PR>>>(n,[&](SizeType i){return Differential<X>::variable(n,deg,z,i)*rad(domain[i]);});
    auto scaled_centre_derivatives = compose(centre_derivatives,scalings);
    auto scaled_derivative_ranges = compose(derivative_ranges,scalings);


    // Make the models
    ValidatedVectorMultivariateTaylorFunctionModel<FLT> tf(m,domain,swp);

    for(SizeType i=0; i!=m; ++i) {
        Differential<FloatBounds<PR>> const& dc = scaled_centre_derivatives[i];
        Differential<FloatBounds<PR>> const& dr = scaled_derivative_ranges[i];
        ValidatedTaylorModel<FLT>& model=tf.model(i);

        Expansion<MultiIndex,FloatDPValue>& expansion=model.expansion();
        FloatError<PR>& error=model.error();
        error=0u;
        expansion.reserve(centre_derivatives[i].expansion().number_of_nonzeros());

        auto riter=dr.begin();
        FloatBounds<PR> coef(pr);

        // Since coefficients are stored in increasing total degree, can first do centre and then ranges
        for(auto centre_iter=dc.begin(); centre_iter!=dc.end() && centre_iter->index().degree()<deg; ++centre_iter) {
            expansion.append(centre_iter->index(),med(centre_iter->coefficient(),error));
        }

        if(not (dr.expansion().empty() or dr.expansion().back().index().degree()<deg)) {
            auto range_iter=dr.begin(); while (range_iter->index().degree()<deg) { ++range_iter; }
            for ( ; range_iter!=dr.end(); ++range_iter) {
                expansion.append(range_iter->index(),med(range_iter->coefficient(),error));
            }
        }

        model.cleanup();
    }
    return tf;
}


// Solve \f$\dt{\phi}(x,t,a)=f(\phi(x,t),t,a)\f$ for x in domx, t in domt, and a in doma, assuming x remains in bndx.
FlowStepModelType
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const ExactIntervalType& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 DegreeType deg,
                 Sweeper<FloatDP> swp)
{
    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension()
                            || f.argument_size()==domx.dimension()+domt.dimension()+doma.dimension());
    const bool is_autonomous = f.argument_size()==domx.dimension()+doma.dimension();

    typedef DoublePrecision PR;
    typedef FloatBounds<PR> X;
    PR pr;

    // Extend time domain from [t:t+h] to [t-h:t+h]
    auto wide_domt = forwards_backwards_time_domain(domt);
    ExactBoxType domxta=product(domx,wide_domt,doma);

    Vector<X> cx(midpoint(domx),pr);
    X t0(domt.lower_bound(),pr);
    Vector<X> ca(midpoint(doma),pr);
    Vector<Differential<X>> cdf = is_autonomous ? f.differential(join(cx,ca),deg) : f.differential(join(cx,t0,ca),deg);
    Vector<Differential<X>> centre_flow_derivatives = is_autonomous ? flow(cdf, cx,ca) : flow(cdf, cx,t0,ca);

    Vector<X> bndx=cast_singleton(bndbx);
    Vector<X> rngx=cast_singleton(domx,pr);
    X rngt=cast_singleton(domt,pr);
    Vector<X> rnga=cast_singleton(doma,pr);
    Vector<Differential<X>> rngdf = is_autonomous ? f.differential(join(bndx,rnga),deg) : f.differential(join(bndx,rngt,rnga),deg);
    Vector<Differential<X>> range_flow_derivatives = is_autonomous ? flow(rngdf, bndx,rnga) : flow(rngdf, bndx,rngt,rnga);

    FlowStepModelType forwards_backwards_taylor_function_model = make_taylor_function_model(domxta, centre_flow_derivatives, range_flow_derivatives, swp);
    domxta[domx.size()]=ExactIntervalType(domt);
    FlowStepModelType forwards_taylor_function_model=restriction(forwards_backwards_taylor_function_model,domxta);
    return forwards_taylor_function_model;
}



FlowStepModelType
TaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const StepSizeType& h, const UpperBoxType& bndx) const
{
    Interval<StepSizeType> domt(0,h);
    ExactBoxType doma;

    return this->flow_step(f, domx,domt,doma, bndx);
}

FlowStepModelType
TaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const Interval<StepSizeType>& rngt, const ExactBoxType& doma, const UpperBoxType& bndx) const
{
    ExactIntervalType domt(rngt);
    FlowStepModelType tphi=Ariadne::series_flow_step(f,domx,domt,doma,bndx,this->order(),this->sweeper());
    return tphi;
}

Pair<StepSizeType,UpperBoxType>
TaylorSeriesIntegrator::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& hmax) const
{
    return configuration().bounder().compute(vf,dx,hmax);
}

Void TaylorSeriesIntegrator::_write(OutputStream& os) const {
    os << "TaylorSeriesIntegrator" << configuration();
}



FlowStepModelType
GradedTaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const StepSizeType& h, const UpperBoxType& bndx) const
{
    Interval<StepSizeType> domt(0,h);
    ExactBoxType doma;

    return this->flow_step(f, domx,domt,doma, bndx);
}

FlowStepModelType
GradedTaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const Interval<StepSizeType>& rngt, const ExactBoxType& doma, const UpperBoxType& bndx) const
{
    ExactIntervalType domt(rngt);
    ExactDouble max_err=this->step_maximum_error();

    DegreeType init_so=this->minimum_spacial_order();
    DegreeType init_to=this->minimum_temporal_order();
    DegreeType max_so=this->maximum_spacial_order();
    DegreeType max_to=this->maximum_temporal_order();

    Vector<ValidatedProcedure> p(f);

    FlowStepModelType tphi=Ariadne::graded_series_flow_step(p,domx,domt,doma,bndx,
        max_err,this->sweeper(), init_so,init_to,max_so,max_to);

    if(tphi.error().raw()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"GradedTaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<domx<<" for time interval "<<domt<<" has error "<<tphi.errors()<<
                      " using spacial order "<<max_so<<" and temporal order "<<max_to<<
                      ", which exceeds maximum single-step error "<<max_err);
    }

    return tphi;
}

Pair<StepSizeType,UpperBoxType>
GradedTaylorSeriesIntegrator::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& hmax) const
{
    return configuration().bounder().compute(vf,dx,hmax);
}

Void GradedTaylorSeriesIntegrator::_write(OutputStream& os) const {
    os << "GradedTaylorSeriesIntegrator" << configuration();
}



template<class X> Void truncate(Differential<X>& x, DegreeType spacial_order_, DegreeType temporal_order_) {
    SizeType n=x.argument_size()-1;
    typename Differential<X>::Iterator write_iter=x.begin();
    typename Differential<X>::ConstIterator read_iter=x.begin();
    while(read_iter!=x.end()) {
        UniformConstReference<MultiIndex> index = read_iter->index();
        if(index[n]>temporal_order_ || index[n]+spacial_order_<index.degree()) {
        } else {
            *write_iter=*read_iter;
            ++write_iter;
        }
        ++read_iter;
    }
    x.expansion().resize(static_cast<SizeType>(write_iter-x.begin()));
}

template<class X> Void truncate(Vector< Differential<X> >& x, DegreeType spacial_order_, DegreeType temporal_order_) {
    for(DegreeType i=0; i!=x.size(); ++i) { truncate(x[i],spacial_order_,temporal_order_); }
}

AffineIntegrator::AffineIntegrator(MaximumError maximum_error_, TemporalOrder temporal_order_)
    : IntegratorBase(/* FIXME */Configuration<TaylorPicardIntegrator>()) { }

AffineIntegrator::AffineIntegrator(MaximumError maximum_error_, SpacialOrder spacial_order_, TemporalOrder temporal_order_)
    : IntegratorBase(/* FIXME */Configuration<TaylorPicardIntegrator>()) { }

Vector<ValidatedDifferential>
AffineIntegrator::flow_derivative(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& dom) const
{
    Vector<ValidatedDifferential> dx=
        ValidatedDifferential::variables(this->_spacial_order+this->_temporal_order,
                                         join(dom,zero));
    dx[dom.size()]=zero;
    Vector<ValidatedDifferential> dphi = dx;
    for(DegreeType i=0; i!=_temporal_order; ++i) {
        dphi = antiderivative(f.evaluate(dphi),dom.size())+dx;
    }
    truncate(dphi,this->_spacial_order,this->_temporal_order);
    return dphi;
}

FlowStepModelType
AffineIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& dom, const StepSizeType& h, const UpperBoxType& bbox) const
{
    Vector<ValidatedNumericType> dmid = Vector<ValidatedNumericType>(midpoint(dom));

    Vector<ValidatedDifferential> mdphi = this->flow_derivative(f,dmid);
    Vector<ValidatedDifferential> bdphi = this->flow_derivative(f,cast_singleton(bbox));

    ARIADNE_WARN("AffineIntegrator may compute overly optimistic error bounds.");

    const SizeType n=dom.size();
    DoublePrecision prec;
    FloatDPError zero_err(prec);

    Vector<FloatDPError> err(n,zero_err);

    Vector<FloatDPError> rad(n+1,zero_err);
    for(SizeType i=0; i!=n; ++i) {
        rad[i] = cast_positive(max(dom[i].upper_bound()-dmid[i].lower(),dmid[i].upper()-dom[i].lower_bound()));
    }
    rad[n] = abs(h);

    for(SizeType i=0; i!=n; ++i) {
        for(Expansion<MultiIndex,ValidatedNumericType>::ConstIterator iter=bdphi[i].begin(); iter!=bdphi[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            if(a[n]==this->_temporal_order && a[n]+this->_spacial_order==a.degree()) {
                UniformConstReference<ValidatedNumericType> rng = iter->coefficient();
                UniformConstReference<ValidatedNumericType> mid = mdphi[i][a];
                ARIADNE_ASSERT(rng.lower().raw()<=mid.lower().raw() && mid.upper().raw()<=rng.upper().raw());
                FloatDPError mag = FloatDPError(max(rng.upper()-mid.lower(),mid.upper()-rng.lower()));
                for(SizeType j=0; j!=n+1u; ++j) { mag *= pow(rad[j],static_cast<DegreeType>(a[j])); }
                err[i] += mag;
            }
        }
    }

    ExactBoxType flow_domain = product(dom,ExactIntervalType(0,h));

    auto function_factory = make_taylor_function_factory();
    FlowStepModelType id = function_factory->create_identity(flow_domain);
    FlowStepModelType res = function_factory->create_zeros(n,flow_domain);
    for(SizeType i=0; i!=n; ++i) {
        ValidatedScalarMultivariateFunctionModelDP res_model = res[i] + mdphi[i].expansion()[MultiIndex::zero(n+1)];
        for(SizeType j=0; j!=mdphi[i].argument_size()-1; ++j) { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*(id[j]-ValidatedNumericType(midpoint(flow_domain[j]))); }
        SizeType j=mdphi[i].argument_size()-1u; { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*id[j]; }
        res_model += FloatDPBounds(-err[i],+err[i]);
        res[i]=res_model;
    }
    return res;
}

FlowStepModelType
AffineIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

Void AffineIntegrator::_write(OutputStream& os) const {
    os << "AffineIntegrator" << configuration();
}



} // namespace Ariadne
