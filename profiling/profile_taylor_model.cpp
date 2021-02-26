/***************************************************************************
 *            profile_taylor_model.cpp
 *
 *  Copyright 2008-21  Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "config.hpp"
#include "function/taylor_model.hpp"
#include "utility/stopwatch.hpp"

using namespace Ariadne;

template<class A1>
struct InplaceUnaryClosure {
    typedef A1 result_type;
    InplaceUnaryClosure(A1&(*f_)(A1&),const A1& a1_) : f(f_), a1(a1_) { };
    A1& operator()() const { r=a1; return f(r); }
    A1&(*f)(A1&); mutable A1 r; const A1& a1;
};

template<class A1, class A2>
struct InplaceBinaryClosure {
    typedef A1 result_type;
    InplaceBinaryClosure(A1&(*f_)(A1&,const A2&),const A1& a1_,const A2& a2_) : f(f_), a1(a1_), a2(a2_) { };
    A1& operator()() const { r=a1; return f(r,a2); }
    A1&(*f)(A1&,const A2&); mutable A1 r; const A1& a1; const A2& a2;
};

template<class R, class A1>
struct UnaryClosure {
    typedef R result_type;
    UnaryClosure(R(*f_)(const A1&),const A1& a1_) : f(f_), a1(a1_) { };
    R operator()() const { return f(a1); }
    R(*f)(const A1&); const A1& a1;
};

template<class R, class A1, class A2>
struct BinaryClosure {
    typedef R result_type;
    BinaryClosure(R(*f_)(const A1&,const A2&),const A1& a1_, const A2& a2_) : f(f_), a1(a1_), a2(a2_) { };
    R operator()() const { return f(a1,a2); }
    R(*f)(const A1&,const A2&); const A1& a1; const A2& a2;
};

template<class F, class A1>
InplaceUnaryClosure<A1>
inplace_bind(const F& f,const A1& a1) {
    return InplaceUnaryClosure<A1>((A1&(*)(A1&))f,a1);
}

template<class F, class A1, class A2>
InplaceBinaryClosure<A1,A2>
inplace_bind(const F& f,const A1& a1,const A2& a2) {
    return InplaceBinaryClosure<A1,A2>((A1&(*)(A1&,const A2&))f,a1,a2);
}

template<class F, class A1>
UnaryClosure<A1,A1>
bind(const F& f,const A1& a1) {
    return UnaryClosure<A1,A1>((A1(*)(const A1&))f,a1);
}

template<class F, class A1, class A2>
BinaryClosure<A1,A1,A2>
bind(const F& f,const A1& a1,const A2& a2) {
    return BinaryClosure<A1,A1,A2>((A1(*)(const A1&,const A2&))f,a1,a2);
}


namespace Ariadne {

typedef TaylorModel<ValidatedTag,FloatDP> TaylorModelType;

}



TaylorModelType copy(const TaylorModelType& x) {
    return x;
}

TaylorModelType& isweep(TaylorModelType& x) {
    x.sweep(); return x;
}

TaylorModelType& bndsadd(TaylorModelType& x, const FloatDPBounds& bnds) {
    return x+=bnds;
}

TaylorModelType& bndsscal(TaylorModelType& x, const FloatDPBounds& bnds) {
    return x*=bnds;
}

TaylorModelType& fscal(TaylorModelType& x, const FloatDPValue& c) {
    return x*=c;
}

TaylorModelType& isum(TaylorModelType& x, const TaylorModelType& y) {
    return x+=y;
}

TaylorModelType sum(const TaylorModelType x1, const TaylorModelType& x2) {
    return x1+x2;
}

TaylorModelType prod(const TaylorModelType x1, const TaylorModelType& x2) {
    return x1*x2;
}

TaylorModelType expo(const TaylorModelType x) {
    return exp(x);
}

TaylorModelType exp_cos(const TaylorModelType x, const TaylorModelType& y) {
    return exp(x)*cos(y);
}

TaylorModelType sigmoid(const TaylorModelType x, const TaylorModelType& y) {
    const FloatDPValue a(10,dp);
    return exp(-x/a);
}

typedef TaylorModelType(*TaylorFunctionPtr)(const Vector<TaylorModelType>&);


template<class T>
Void profile(Nat ntries, const char* name, const T& run)
{
    typename T::result_type res=run();
    //std::cerr<< "\n" << name << "(" << args << ")=\n  " << res << "\n\n";

    StopWatch sw;

    for(Nat i=0; i!=ntries; ++i) {
        res=run();
    }

    sw.click();

    TimeUnit total_time = sw.elapsed();
    double average_time_in_microseconds = 1000000*total_time/(double)ntries;
    FloatDPError error = res.error();
    SizeType size = res.number_of_nonzeros();

    std::cout << std::setw(20) << std::left << name << std::right
              << std::setw(10) << std::fixed << std::setprecision(2) << average_time_in_microseconds << " "
              << std::setw(12) << std::scientific << std::setprecision(4) << error.raw().get_d()
              << std::setw(8) << size
              << std::endl;
}



Int main(Int argc, const char* argv[]) {
    Nat nbase=10;
    if(argc>1) { nbase=atoi(argv[1]); }

    Vector<FloatDPValue> c({1.0_x,2.0_x},dp);
    Int i;

    TrivialSweeper<FloatDP> trivial_sweeper(dp);
    Sweeper<FloatDP> threshold_sweeper=ThresholdSweeper<FloatDP>(dp,1e-5);

    // Use for in clean()
    Vector<TaylorModelType> v(2,TaylorModelType(2,threshold_sweeper));
    v[0]=TaylorModelType::coordinate(2,0,threshold_sweeper);
    v[1]=TaylorModelType::constant(2,1.0_x,threshold_sweeper);

    // Use in clean()
    TaylorModelType w(3,trivial_sweeper);
    for(MultiIndex a(3); a.degree()<=9; ++a) {
        if(i%7<3) { w.expansion().append(a,cast_exact(1/(1.0+i*i*i*i*i))); }
        else if(i%7<4) { w.expansion().append(a,cast_exact(1/(1.0+i))); }
        ++i;
    }
    w.set_sweeper(threshold_sweeper);

    // Use in arithmetical operations
    TaylorModelType x(3,trivial_sweeper);
    TaylorModelType y(3,trivial_sweeper);
    i=0;
    for(MultiIndex a(3); a.degree()<=7; ++a) {
        if(i%7<4) { x.expansion().append(a,cast_exact(1/(1.0+i))); }
    }
    i=0;
    for(MultiIndex a(3); a.degree()<=5; ++a) {
        if(i%7<4) { x.expansion().append(a,cast_exact(1/(1.0+i))); }
        if(i%3<2) { y.expansion().append(a,cast_exact(1/(2.0+i))); }
        ++i;
    }
    x.cleanup();
    y.cleanup();

    Expansion<MultiIndex,FloatDPValue> zc({ {{0,0,0},1.0_x}, {{1,0,0},0.5_x}, {{0,1,0},-0.25_x}, {{0,0,1},0.625_x}, {{1,1,0},0.375_x} },dp);
    FloatDPError ze(dp);
    TaylorModelType z(zc,ze,threshold_sweeper);
    z.cleanup();

    FloatDPBounds bnds(Rational(0.33_dec),0.49_dec,dp);
    FloatDPValue cnst(0.41015625_x,dp);

    std::cout << std::setw(20) << std::left << "name" << std::right
              << std::setw(11) << "time(us)"
              << std::setw(12) << "error"
              << std::setw(8) << "size"
              << std::endl;

    profile(nbase*10000,"isweep-02",inplace_bind(&isweep,w));
    profile(nbase*10000,"copy-02",bind(&copy,w));
    profile(nbase*10000,"iadd-noinsert-02",inplace_bind(&bndsadd,x,bnds));
    profile(nbase*10000,"iadd-insert-02",inplace_bind(&bndsadd,x,bnds));
    profile(nbase*10000,"iscal-02",inplace_bind(&bndsscal,x,bnds));
    profile(nbase*10000,"fscal-02",inplace_bind(&fscal,x,cnst));
    profile(nbase*1000,"isum-02",inplace_bind(&isum,x,y));
    profile(nbase*1000,"sum-02",bind(&sum,x,y));
    profile(nbase*100,"prod-02",bind(&prod,x,y));
    profile(nbase*100,"exp-02",bind(&expo,z));
    profile(nbase*100,"exp_cos-01",bind(exp_cos,v[0],v[1]));
    profile(nbase*100,"sigmoid-01",bind(sigmoid,v[0]));
}
