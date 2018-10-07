/***************************************************************************
 *            concepts.hpp
 *
 *  Copyright  2018  Pieter Collins
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

/*! \file concepts.hpp
 *  \brief Function concept classes.
 */

#ifndef ARIADNE_CONCEPTS_HPP
#define ARIADNE_CONCEPTS_HPP

#include "../utility/metaprogramming.hpp"

namespace Ariadne {

template<class X> class ConcreteConcept {
  public:
    typedef typename X::GenericType Y;
    typedef typename X::PropertiesType PR;

    void usage() {
        { PR pr = x.properties(); }
        { Y y = x.generic(); }
        { x.operator Y(); }
        { X x(pr); }
        { X x(y,pr); }
        { X x=create(y,pr); }
        { x_nc=y; }
    }
  private:
    X x_nc;
    X const x;
    Y const y;
    PR const pr;
};

class ConcreteArchetype {
  private:
    typedef ConcreteArchetype X;
    struct Y { };
    struct PR { };
    friend X create(Y const&, PR);
  public:
    typedef Y GenericType;
    typedef PR PropertiesType;
    PR properties() const;
    Y generic() const;
    operator Y() const;
    ConcreteArchetype(PR);
    ConcreteArchetype(Y const&, PR);
    X& operator=(Y const&);
};
ConcreteArchetype::X create(ConcreteArchetype::Y const&, ConcreteArchetype::PR);

template<class R> class FieldConcept {
  private:
    R r;
  public:
    void usage() {
        R rr(r);
        r=neg(r); 
        r=rec(r);
        r=add(r,r); 
        r=sub(r,r);
        r=mul(r,r);
        r=div(r,r);
        
        r=+r; r=-r; 
        r=r+r; r=r-r; r=r*r; r=r/r;
        r+=r; r-=r; r*=r; r/=r;
    }
};

template<class R> class ElementaryConcept {
  private:
    R r;
  public:
    void usage() {
        r=sqrt(r); 
        r=exp(r);
        r=log(r);
        r=sin(r);
        r=cos(r);
        r=tan(r);
        r=atan(r);
    }
};

template<class R> class CompleteFieldConcept : public FieldConcept<R>, public ElementaryConcept<R> {
};

class Rational;

template<class R> class RealConcept
    : public CompleteFieldConcept<R>
{
  private:
    Rational* q;

    void usage() {
        R r(*q);
    }
};

template<class A> class AlgebraConcept
    : FieldConcept<A>
{
  public:
    typedef typename A::NumericType X;

    void usage() {
        a_nc = x;
        a = add(x,a);
        a = sub(x,a);
        a = mul(x,a);
        a = div(x,a);
        a = add(a,x);
        a = sub(a,x);
        a = mul(a,x);
        a = div(a,x);
        
        a=x+a; a=x-a; a=x*a; a=x/a;
        a=a+x; a=a-x; a=a*x; a=a/x;
        a+=x; a-=x; a*=x; a/=x;        
    }
  private:
    A a_nc;
    A const a;
    X const x;
};

class RealArchetype {
  private:
    typedef RealArchetype Self;
  public:
    RealArchetype(Self const&);
    RealArchetype(Rational const&);
    friend Self neg(Self const&);
    friend Self add(Self const&, Self const&);
    friend Self sub(Self const&, Self const&);
    friend Self mul(Self const&, Self const&);
    friend Self div(Self const&, Self const&);
    friend Self rec(Self const&);
};

template<class F> class FunctionConcept {
    typedef typename F::Paradigm P;
    template<class A> using Argument = typename F::template Argument<A>;
    template<class A> using Result = typename F::template Result<A>;

    F f; SizeType k;
    Argument<Number<P> ay; Argument<FloatDPBounds> ax; 
    Result<Number<P>> ry; Result<FloatDPBounds> rx; 
  public:
    bool check() {
        f(ay);
        ry=evaluate(f,ay); 
        ry=unchecked_evaluate(f,ay);
        f=partial_evaluate(f,k,y);
        
        if constexpr (IsSame<P,ValidatedTag>::value) {
            f(ax);
            rx=evaluate(f,ax);
            rx=unchecked_evaluate(f,ax);
            f=partial_evaluate(f,k,x);
        }

        //a=antiderivative(a,k); a=antiderivative(a,k,x);

        return true;
    }
};
    
    

template<class A> class CompleteAlgebraConcept 
    : public AlgebraConcept<A>, public ElementaryConcept<A>
{
    bool usage() {
        norm(a);
        return true;
    }
};

template<class A> class ConcreteAlgebraConcept : public CompleteAlgebraConcept {
    typedef typename A::PrecisionType PR;
    typedef typename A::GenericType B;
    
    typedef typename A::NumericType X;
    typedef typename X::GenericType Y;

    static_assert(IsSame<typename X::GenericType,typename A::GenericNumericType>::value);
    
    B b; Y y;
    
  public:    
    bool usage() {
        a=b; a=y;
        
        a=a+b; a=a-b; a=a*b; a=a/b;
        a=b+a; a=b-a; a=b*a; a=b/a;
        
        a=a+y; a=a-y; a=a*y; a=a/y;
        a=y+a; a=y-a; a=y*a; a=y/a;
        
        a+=b; a-=b; a*=b; a/=b;
        a+=y; a-=y; a*=y; a/=y;
        
        a=add(a,b); a=sub(a,b); a=mul(a,b); a=div(a,b);
        a=add(b,a); a=sub(b,a); a=mul(b,a); a=div(b,a);
        
        a=add(a,y); a=sub(a,y); a=mul(a,y); a=div(a,y);
        a=add(y,a); a=sub(y,a); a=mul(y,a); a=div(y,a);

        return true;
    }
};

} // namespace Ariadne 

#endif // ARIADNE_CONCEPTS_HPP
