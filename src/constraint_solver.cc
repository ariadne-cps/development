/***************************************************************************
 *            constraint_solver.cc
 *
 *  Copyright 2000  Pieter Collins
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

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

#include "macros.h"
#include "tuple.h"
#include "tribool.h"
#include "numeric.h"
#include "vector.h"
#include "box.h"
#include "grid_set.h"
#include "polynomial.h"
#include "function.h"
#include "formula.h"
#include "procedure.h"
#include "constraint.h"
#include "nonlinear_programming.h"
#include "function_mixin.h"
#include "taylor_function.h"

#include "constraint_solver.h"
#include "solver.h"

namespace Ariadne {


static const double error =  1e-2;

enum Sign { NEGATIVE=-1, ZERO=0, POSITIVE=+1 };

Sign sign(const Float& x) {
    if(x>0) { return NEGATIVE; }
    else if(x<0) {  return POSITIVE; }
    else { return ZERO; }
}

Sign sign(const Interval& ivl) {
    if(ivl.lower()>0) { return NEGATIVE; }
    else if(ivl.upper()<0) {  return POSITIVE; }
    else { return ZERO; }
}

typedef Vector<Float> FloatVector;
typedef Matrix<Float> FloatMatrix;
typedef boost::numeric::ublas::vector_range<FloatVector> FloatVectorRange;

typedef Vector<Interval> IntervalVector;
typedef Matrix<Interval> IntervalMatrix;
typedef boost::numeric::ublas::vector_range<IntervalVector> IntervalVectorRange;


std::ostream& operator<<(std::ostream& os, const NonlinearConstraint& c) {
    static const Float inf = Ariadne::inf<Float>();
    if(c.bounds().lower()==c.bounds().upper()) { return os << c.function() << "==" << c.bounds().upper(); }
    if(c.bounds().upper()==inf) { return os << c.bounds().lower() << "<=" << c.function(); }
    if(c.bounds().lower()==-inf) { return os << c.function() << "<=" << c.bounds().upper(); }
    return os << c.bounds().lower() << "<=" << c.function() << "<=" << c.bounds().upper();
}



Pair<Tribool,Point> ConstraintSolver::feasible(const Box& domain, const List<NonlinearConstraint>& constraints) const
{
    if(constraints.empty()) { return make_pair(!domain.empty(),domain.centre()); }

    RealVectorFunction function(constraints.size());
    Box bounds(constraints.size());

    for(uint i=0; i!=constraints.size(); ++i) {
        function.set(i,constraints[i].function());
        bounds[i]=constraints[i].bounds();
    }
    return this->feasible(domain,function,bounds);
}


Pair<Tribool,Point> ConstraintSolver::feasible(const Box& domain, const IntervalVectorFunction& function, const Box& codomain) const
{

    static const double XSIGMA=0.125;
    static const double TERR=-1.0/(1<<10);
    static const Float inf = Ariadne::inf<Float>();

    ARIADNE_LOG(4,"domain="<<domain<<"\nfunction="<<function<<"\ncodomain="<<codomain<<"\n");

    // Make codomain bounded
    IntervalVector bounds=codomain;
    IntervalVector image=function(domain);
    ARIADNE_LOG(4,"image="<<image<<"\n");
    for(uint i=0; i!=image.size(); ++i) {
        if(disjoint(image[i],codomain[i])) { ARIADNE_LOG(4,"  Proved disjointness using direct evaluation\n");  return make_pair(false,Point()); }
        else bounds[i]=intersection(codomain[i],image[i]);
    }


    const uint m=domain.size(); // The total number of variables
    const uint n=codomain.size(); // The total number of nontrivial constraints
    const uint l=(m+n)*2; // The total number of lagrange multipliers

    Point point(m); // The point in the domain which is the current test point
    Float violation; // An upper bound on amount by which the constraints are violated by the test point
    Point multipliers(l); // The lagrange multipliers for the constraints
    Point slack(l); // The slack between the test point and the violated constraints

    Float& t=violation; Point& x=multipliers; Point& y=point; Point& z=slack; // Aliases for the main quantities used
    const Box& d=domain; const IntervalVectorFunction& fn=function; const Box& c=codomain; // Aliases for the main quantities used
    VectorTaylorFunction tfn(d,fn);

    point=midpoint(d);
    for(uint k=0; k!=l; ++k) { multipliers[k]=1.0/l; }

    NonlinearInteriorPointOptimiser optimiser;
    optimiser.compute_tz(domain,function,bounds,point,violation,slack);

    ARIADNE_LOG(4,"d="<<d<<", f="<<fn<<", c="<<c<<"\n");


    // TODO: Don't use fixed number of steps
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"    t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
        optimiser.feasibility_step(d,fn,c,x,y,z,t);
        if(t>=TERR) {
            ARIADNE_LOG(4,"t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
            if(this->check_feasibility(domain,function,codomain,point)) { return make_pair(true,point); }
            else { ARIADNE_LOG(2,"f(y)="<<fn(IntervalVector(y))<<"\n"); return make_pair(indeterminate,point); }
        }
    }
    ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        Box subdomain=domain;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ScalarTaylorFunction txg=ScalarTaylorFunction::constant(d,0);
        Interval cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            txg = txg - Interval(x[j]-x[n+j])*tfn[j];
            cnst += (c[j].upper()*x[j]-c[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            txg = txg - (x[2*n+i]-x[2*n+m+i])*ScalarTaylorFunction::coordinate(d,i);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        txg = cnst + txg;

        ARIADNE_LOG(4,"    txg="<<txg<<"\n");

        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        this->hull_reduce(subdomain,txg,Interval(0,inf));
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        if(subdomain.empty()) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return make_pair(false,Point());
        }

        for(uint i=0; i!=m; ++i) {
            this->box_reduce(subdomain,txg,Interval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<subdomain<<"\n");
            if(subdomain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return make_pair(false,Point()); }
        }
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");

        //Pair<Box,Box> sd=solver.split(List<NonlinearConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=d.split();
        Point nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point ny = midpoint(sd.first);
        tribool result=this->feasible(sd.first, fn, c).first;
        nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        result = result || this->feasible(sd.second, fn, c).first;
        return make_pair(result,Point());
    }

    return make_pair(indeterminate,Point());
}


void ConstraintSolver::reduce(Box& domain, const IntervalVectorFunction& function, const Box& codomain) const
{
    const double MINIMUM_REDUCTION = 0.75;
    ARIADNE_ASSERT(function.argument_size()==domain.size());
    ARIADNE_ASSERT(function.result_size()==codomain.size());

    if(domain.empty()) { return; }

    Float domain_magnitude=0.0;
    for(uint j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width();
    }
    Float old_domain_magnitude=domain_magnitude;

    do {
        this->hull_reduce(domain,function,codomain);
        if(domain.empty()) { return; }

        for(uint i=0; i!=codomain.size(); ++i) {
            for(uint j=0; j!=domain.size(); ++j) {
                this->box_reduce(domain,function[i],codomain[i],j);
            }
        }
        if(domain.empty()) { return; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0.0;
        for(uint j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude/old_domain_magnitude <= MINIMUM_REDUCTION);
}

void ConstraintSolver::reduce(Box& domain, const List<NonlinearConstraint>& constraints) const
{
    const double MINIMUM_REDUCTION = 0.75;

    if(domain.empty()) { return; }

    Float domain_magnitude=0.0;
    for(uint j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width();
    }
    Float old_domain_magnitude=domain_magnitude;

    do {
        for(uint i=0; i!=constraints.size(); ++i) {
            this->hull_reduce(domain,constraints[i].function(),constraints[i].bounds());
        }

        if(domain.empty()) { return; }

        for(uint i=0; i!=constraints.size(); ++i) {
            for(uint j=0; j!=domain.size(); ++j) {
                this->box_reduce(domain,constraints[i].function(),constraints[i].bounds(),j);
            }
        }
        if(domain.empty()) { return; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0.0;
        for(uint j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude/old_domain_magnitude <= MINIMUM_REDUCTION);
}


void ConstraintSolver::hull_reduce(Box& domain, const IntervalProcedure& procedure, const Interval& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, IntervalProcedure procedure, Interval bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
}

void ConstraintSolver::hull_reduce(Box& domain, const Vector<IntervalProcedure>& procedure, const IntervalVector& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, Vector<IntervalProcedure> procedure, Box bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
}

void ConstraintSolver::hull_reduce(Box& domain, const IntervalScalarFunctionInterface& function, const Interval& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, IntervalScalarFunction function, Interval bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Formula<Interval> formula=function.evaluate(Formula<Interval>::identity(function.argument_size()));
    Procedure<Interval> procedure(formula);
    this->hull_reduce(domain,procedure,bounds);
}

void ConstraintSolver::hull_reduce(Box& domain, const IntervalVectorFunctionInterface& function, const IntervalVector& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, IntervalScalarFunction function, Interval bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Vector< Formula<Interval> > formula=function.evaluate(Formula<Interval>::identity(function.argument_size()));
    Vector< Procedure<Interval> > procedure(formula);
    this->hull_reduce(domain,procedure,bounds);
}

void ConstraintSolver::monotone_reduce(Box& domain, const IntervalScalarFunctionInterface& function, const Interval& bounds, uint variable) const
{
    IntervalScalarFunction derivative=function.derivative(variable);

    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<", derivative="<<derivative<<"\n");

    Float midpoint;
    Interval lower=domain[variable];
    Interval upper=domain[variable];
    Vector<Interval> slice=domain;
    Vector<Interval> subdomain=domain;

    static const int MAX_STEPS=3;
    const Float size = lower.radius() / (1<<MAX_STEPS);
    do {
        if(lower.radius()>size) {
            midpoint=lower.midpoint();
            slice[variable]=midpoint;
            Interval new_lower=midpoint+(bounds-function(slice)).lower()/derivative(subdomain);
            if(new_lower.upper()<lower.lower()) { lower=lower.lower(); }
            else { lower=intersection(lower,new_lower); }
        }
        if(upper.radius()>size) {
            midpoint=upper.midpoint();
            slice[variable]=midpoint;
            Interval new_upper=midpoint+(bounds-function(slice)).upper()/derivative(subdomain);
            if(new_upper.lower()>upper.upper()) { upper=upper.upper(); }
            else { upper=intersection(upper,new_upper); }
        }
        subdomain[variable]=Interval(lower.lower(),upper.upper());
    } while(lower.radius()>size && upper.radius()>size);
    domain=subdomain;
}





void ConstraintSolver::box_reduce(Box& domain, const IntervalScalarFunctionInterface& function, const Interval& bounds, uint variable) const
{
    ARIADNE_LOG(2,"ConstraintSolver::box_reduce(Box domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<"\n");

    Interval interval=domain[variable];
    Float l=interval.lower();
    Float u=interval.upper();
    Interval subinterval;
    Interval new_interval(+inf<Float>(),-inf<Float>());
    Vector<Interval> slice=domain;

    static const uint MAX_STEPS=(1<<3);
    const uint n=MAX_STEPS;
    for(uint i=0; i!=n; ++i) {
        subinterval=Interval((l*(n-i)+u*i)/n,(l*(n-i-1)+u*(i+1))/n);
        slice[variable]=subinterval;
        Interval slice_image=function(slice);
        if(!intersection(slice_image,bounds).empty()) {
            new_interval.set_lower(subinterval.lower());
            break;
        }
    }

    for(uint i=0; i!=n; ++i) {
        subinterval=Interval((u*(n-i-1)+l*(i+1))/n,(u*(n-i)+l*i)/n);
        slice[variable]=subinterval;
        Interval slice_image=function(slice);
        if(!intersection(slice_image,bounds).empty()) {
            new_interval.set_upper(subinterval.upper());
            break;
        }
    }
    domain[variable]=new_interval;
}



namespace {

void compute_monotonicity(Box& domain, const NonlinearConstraint& constraint) {
/*
    static const uint n = domain.size();

    // Compute monotone formulae
    boost::array<Sign,0> monotonicity(n);
    Vector<Interval> grad=constraint.function().gradient(domain);
    for(uint j=0; j!=n; ++j) {
        monotonicity[j]=sign(grad[j]);
    }
*/
}
} // namespace


Pair<Box,Box> ConstraintSolver::split(const Box& d, const IntervalVectorFunction& f, const IntervalVector& c) const
{
    return d.split();
}


tribool ConstraintSolver::check_feasibility(const Box& d, const IntervalVectorFunction& f, const Box& c, const Point& y) const
{
    for(uint i=0; i!=y.size(); ++i) {
        if(y[i]<d[i].lower() || y[i]>d[i].upper()) { return false; }
    }

    Vector<Interval> fy=f(Vector<Interval>(y));
    ARIADNE_LOG(4,"d="<<d<<" f="<<f<<", c="<<c<<"\n  y="<<y<<", f(y)="<<fy<<"\n");
    tribool result=true;
    for(uint j=0; j!=fy.size(); ++j) {
        if(fy[j].lower()>c[j].upper() || fy[j].upper()<c[j].lower()) { return false; }
        if(fy[j].upper()>=c[j].upper() || fy[j].lower()<=c[j].lower()) { result=indeterminate; }
    }
    return result;
}






} // namespace Ariadne
