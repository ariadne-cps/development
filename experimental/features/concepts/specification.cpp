#include <type_traits>
#include <functional>
#include <initializer_list>
#include <vector>
#include <set>
#include <map>


using Void = void;
using SizeType = std::size_t;

using True = std::true_type;

template<class T> using InitializerList = std::initializer_list<T>;
template<class T> using Array = std::vector<T>;
template<class T> using List = std::vector<T>;
template<class K, class V> using Map = std::map<K,V>;
template<class SIG> using Function = std::function<SIG>;

using std::declval;
template<class P, class V=Void> using EnableIf = typename std::enable_if<P::value,V>::type;

template<class I, class X> class Expansion : public Map<I,X> {
    typedef I IndexType;
    typedef X CoefficientType;

    typedef typename Map<I,X>::iterator Iterator;
    typedef typename Map<I,X>::const_iterator ConstIterator;
    
    X& operator[](IndexType i) { return this->Map<I,X>::operator[](i); }
    X const& operator[](IndexType i) const { return const_cast<Map<I,X>*>(this)->Map<I,X>::operator[](i); }

    Void append(IndexType, CoefficientType);
};


struct ExactInformation { };
struct RoundedInformation { };

template<class X> class Effective;
template<class X> class Validated;
template<class X> class Approximate;


class Effort {
    Effort(int n) = delete;
  public:
    Effort(uint m);
};


class Sierpinskian;

template<> class Effective<Sierpinskian> {
  public:
    Validated<Sierpinskian> check(Effort) const;
    friend Validated<Sierpinskian> check(Effective<Sierpinskian> s, Effort eff);
};
template<> class Validated<Sierpinskian> {
    friend bool definitely(Validated<Sierpinskian>);
    friend Validated<Sierpinskian> check(Effective<Sierpinskian> s, Effort eff) { return s.check(eff); }
};


template<class X> using NegType = decltype(neg(declval<X>()));
template<class X, class Y=X> using AddType = decltype(add(declval<X>(),declval<Y>()));
template<class X, class Y=X> using SubType = decltype(sub(declval<X>(),declval<Y>()));
template<class X, class Y=X> using MulType = decltype(mul(declval<X>(),declval<Y>()));
template<class X, class Y=X> using DivType = decltype(div(declval<X>(),declval<Y>()));

template<class X, class Y=X> struct ArithmeticOperations {
    friend X add(X,Y); friend X sub(X,Y); friend X mul(X,Y); friend X div(X,Y);
    friend X add(Y,X); friend X sub(Y,X); friend X mul(Y,X); friend X div(Y,X);
};
template<class X> struct ArithmeticOperations<X> {
    friend X neg(X); friend X hlf(X); friend X add(X,X); friend X sub(X,X); friend X mul(X,X); friend X div(X,X);
};

template<class X> class Vector;

template<class A, class X> class AlgebraOperations {
    friend A add(A,A); friend A sub(A,A); friend A mul(A,A);
    friend A add(X,A); friend A sub(X,A); friend A mul(X,A);
    friend A add(A,X); friend A sub(A,X); friend A mul(A,X); friend A div(A,X);    
};

template<class X> class VectorTraits;
template<class X> using ScalarType = typename VectorTraits<X>::ScalarType;
template<class X> struct VectorTraits<Vector<X>> { typedef X ScalarType; };

template<class VA, class X, template<class>class P> class AlgebraOperations<VA,P<Vector<X>>> {
    using SA=ScalarType<VA>; 
    using SX=P<X>;
    using VX=Vector<P<X>>;
  public:
    SA operator[](SizeType i) const;
    friend VA add(VA,VA); friend VA sub(VA,VA); friend VA mul(SA,VA); friend VA mul(VA,SA);
    friend VA add(VA,VX); friend VA sub(VA,VX); friend VA mul(VA,SX); friend VA div(VA,SX); 
    friend VA add(VX,VA); friend VA sub(VX,VA); friend VA mul(SX,VA); 
};

class Integer; class Dyadic; class Rational;
template<> struct ArithmeticOperations<Integer> {
    using Z=Integer; using W=Dyadic; using Q=Rational;
    friend Z neg(Z); friend W hlf(Z); friend Z add(Z,Z); friend Z sub(Z,Z); friend Z mul(Z,Z); friend Q div(Z,Z);
};
template<> struct ArithmeticOperations<Dyadic> {
    using W=Dyadic; using Q=Rational;
    friend W neg(W); friend W hlf(W); friend W add(W,W); friend W sub(W,W); friend W mul(W,W); friend Q div(W,W);
};
template<> struct ArithmeticOperations<Rational> {
    using Q=Rational;
    friend Q neg(Q); friend Q hlf(Q); friend Q add(Q,Q); friend Q sub(Q,Q); friend Q mul(Q,Q); friend Q div(Q,Q);
};




class Dyadic : ArithmeticOperations<Dyadic> {
  public:
    typedef ExactInformation InformationTag;
    Dyadic(int, uint=0u); Dyadic(double)=delete;
};
    
class Rational : ArithmeticOperations<Rational> {
    Rational(double) = delete;
  public:
    Rational(int, int=1);
    Rational(Dyadic);
};


template<class X> class Error;
template<class X> class Value;
template<class X> class Bounds;
template<class X> class Approximation;

// FIXME: Check for type
template<> class Bounds<Dyadic> {
    using X = Dyadic;
  protected:
    X _l, _u;
  public:
    explicit Bounds(X l, X u) : _l(l), _u(u) { }
    friend Bounds<X> add(Bounds<X> x1, Bounds<X> x2) {
        return Bounds<X>(add(x1._l,x2._l),add(x1._u,x2._u)); }
    friend Bounds<X> sub(Bounds<X> x1, Bounds<X> x2) {
        return Bounds<X>(sub(x1._l,x2._u),sub(x1._u,x2._l)); }
    friend Bounds<X> mul(Bounds<X> x1, Bounds<X> x2);
};


/*
template<> class Approximation<Dyadic> {
    using X = Dyadic;
  protected:
    X _a;
  public:
    explicit Approximation(X a) : _a(a) { }
    friend Approximation<X> add(Approximation<X> x1, Approximation<X> x2) {
        return Approximation<X>(add(x1._a,x2._a)); }
    friend Approximation<X> sub(Approximation<X> x1, Approximation<X> x2) {
        return Approximation<X>(sub(x1._a,x2._a)); }
    friend Approximation<X> mul(Approximation<X> x1, Approximation<X> x2) {
        return Approximation<X>(mul(x1._a,x2._a)); }
};
*/



class DoublePrecision;
class MultiplePrecision;


class Real : ArithmeticOperations<Real> {
  public:
    typedef Real GenericType;
    Real(Dyadic);
    Real(Rational);
};

template<> class Effective<Real> : ArithmeticOperations<Effective<Real>> {
  public:
    typedef Effective<Real> GenericType;
    Effective<Real>(Rational);
};
    
template<> class Validated<Real> : ArithmeticOperations<Validated<Real>> {
  public:
    typedef Validated<Real> GenericType;
    Validated<Real>(Effective<Real>, Effort);
    Bounds<Dyadic> get() const;
};
    
template<> class Approximate<Real> : ArithmeticOperations<Approximate<Real>> {
  public:
    Approximate<Real>(Effective<Real>, Effort);
    Approximate<Real>(Validated<Real>);
    Approximation<Dyadic> get() const;
};

using EffectiveReal = Effective<Real>;
using ValidatedReal = Validated<Real>;
using ApproximateReal = Approximate<Real>;

enum class RoundingMode : char { DOWN=-1, NEAR=0, UP=+1 };
static const RoundingMode down=RoundingMode::DOWN;
static const RoundingMode near=RoundingMode::NEAR;
static const RoundingMode up=RoundingMode::UP;


template<class PR> class Float {
    using RND = RoundingMode;
    using FLT = Float<PR>;
  public:  
    typedef RoundedInformation InformationTag;
    
    typedef PR PrecisionType;
    typedef RND RoundingModeType;
    Float(Dyadic, PR);
    Float(Rational, RND, PR);
    friend FLT neg(FLT);
    friend FLT hlf(FLT);
    friend FLT abs(FLT);
    friend FLT min(FLT,FLT);
    friend FLT max(FLT,FLT);
    friend FLT add(RND, FLT, FLT);
    friend FLT sub(RND, FLT, FLT);
    friend FLT mul(RND, FLT, FLT);
    friend FLT div(RND, FLT, FLT);
};

template<class PR> class Float;
template<class PR> using FloatBounds = Bounds<Float<PR>>;

template<class T> using PrecisionType = typename T::PrecisionType;

class DoublePrecision { 
};

class MultiplePrecision { 
    MultiplePrecision(uint m);
};


using DoublePrecisionFloat = Float<DoublePrecision>;
using MultiplePrecisionFloat = Float<MultiplePrecision>;

template<class X> class Operations {
  public:
    static Bounds<X> _mul(Bounds<X>, Bounds<X>);
    static Bounds<X> _div(Bounds<X>, Bounds<X>);
};

template<class X> class PositiveUpperBound {
  protected:
    X _e;
  public:
    friend PositiveUpperBound<X> add(PositiveUpperBound<X> x1, PositiveUpperBound<X> x2) {
        return PositiveUpperBound<X>(add(up,x1._e,x2._e)); }
    friend PositiveUpperBound<X> mul(PositiveUpperBound<X> x1, PositiveUpperBound<X> x2) {
        return PositiveUpperBound<X>(mul(up,x1._e,x2._e)); }
    friend PositiveUpperBound<X> max(PositiveUpperBound<X> x1, PositiveUpperBound<X> x2) {
        return PositiveUpperBound<X>(max(x1._e,x2._e)); }
    friend PositiveUpperBound<X> min(PositiveUpperBound<X> x1, PositiveUpperBound<X> x2) {
        return PositiveUpperBound<X>(min(x1._e,x2._e)); }
};

template<class X> using ErrorBound = PositiveUpperBound<X>;

template<class X> class Value {
    friend class Bounds<X>;
    friend class Approximation<X>;
  protected:
    X _v;
  public:
    Value(Dyadic w, PrecisionType<X> pr) : _v(w,pr) { };
    friend Value<X> neg(Value<X> vx) {
        return Value<X>(neg(vx._v)); }
    template<class XE> friend Value<X> add(Value<X> x1, Value<X> x2, Error<XE>& e) {
        Value<X> r=add(near,x1,x2); Error<X> ne=hlf(sub(up,add(up,x1,x2),add(down,x1,x2))); e=add(up,e,ne); return r; }
    template<class XE> friend Value<X> sub(Value<X> x1, Value<X> x2, Error<XE>& e) {
        Value<X> r=sub(near,x1,x2); Error<X> ne=hlf(sub(up,sub(up,x1,x2),sub(down,x1,x2))); e=add(up,e,ne); return r; }
    template<class XE> friend Value<X> mul(Value<X> x1, Value<X> x2, Error<XE>& e) {
        Value<X> r=mul(near,x1,x2); Error<X> ne=hlf(sub(up,mul(up,x1,x2),mul(down,x1,x2))); e=add(up,e,ne); return r; }
};    
    
template<class X> class Bounds {
    friend class Approximation<X>;
  protected:
    X _l, _u;
  public:
    explicit Bounds(X l, X u) : _l(l), _u(u) { }
    Bounds(Value<X> vx) : _l(vx._v), _u(vx._v) { }
    Bounds(Rational q, PrecisionType<X> pr) : _l(q,down,pr), _u(q,up,pr) { };
    Bounds(Validated<Real>, PrecisionType<X>);
    operator Validated<Real> ()  const;
    friend Bounds<X> neg(Bounds<X> vx) {
        return Bounds<X>(neg(vx._u),neg(vx._l)); }
    friend Bounds<X> hlf(Bounds<X> vx) {
        return Bounds<X>(hlf(vx._l),hlf(vx._u)); }
    friend Bounds<X> add(Bounds<X> x1, Bounds<X> x2) {
        return Bounds<X>(add(down,x1._l,x2._l),add(up,x1._u,x2._u)); }
    friend Bounds<X> sub(Bounds<X> x1, Bounds<X> x2) {
        return Bounds<X>(sub(down,x1._l,x2._u),sub(up,x1._u,x2._l)); }
    friend Bounds<X> mul(Bounds<X> x1, Bounds<X> x2) {
        return Operations<X>::_mul(x1,x2); }
    friend Bounds<X> div(Bounds<X> x1, Bounds<X> x2) {
        return Operations<X>::_div(x1,x2); }
};

template<class Y, class PR> Bounds<Float<PR>> float_bounds(Y y, PR pr) { return Bounds<Float<PR>>(y,pr); }



template<class X> class Approximation {
    static const bool _rnd = std::is_same<typename X::InformationTag, RoundedInformation>::value;
    static X _med(X l, X u) { if constexpr (_rnd) { return hlf(add(near,l,u)); } else { return hlf(add(l,u)); } }
  protected:
    X _a;
  public:
    explicit Approximation(X a) : _a(a) { }
//    Approximation(Approximate<Real>, PrecisionType<X>);
    Approximation(Value<X> vx) : _a(vx._v) { }
    Approximation(Bounds<X> vx) : _a(_med(vx._l,vx._u)) { }
    friend Approximation<X> neg(Approximation<X> vx) {
        return Approximation<X>(neg(vx._a)); }
    friend Approximation<X> hlf(Approximation<X> vx) {
        return Approximation<X>(hlf(vx._a)); }
    friend Approximation<X> add(Approximation<X> x1, Approximation<X> x2) {
        if constexpr (_rnd) { return Approximation<X>(add(near,x1._a,x2._a)); }
        else { return Approximation<X>(add(x1._a,x2._a)); } }
    friend Approximation<X> sub(Approximation<X> x1, Approximation<X> x2) {
        if constexpr (_rnd) { return Approximation<X>(sub(near,x1._a,x2._a)); }
        else { return Approximation<X>(sub(x1._a,x2._a)); } }
    friend Approximation<X> mul(Approximation<X> x1, Approximation<X> x2) {
        if constexpr (_rnd) { return Approximation<X>(mul(near,x1._a,x2._a)); }
        else { return Approximation<X>(mul(x1._a,x2._a)); } }
    friend Approximation<X> div(Approximation<X> x1, Approximation<X> x2) {
        if constexpr (_rnd) { return Approximation<X>(div(near,x1._a,x2._a)); }
        else { return Approximation<X>(div(x1._a,x2._a)); } }
};


typedef unsigned short Index;
typedef Array<Index> MultiIndex;


template<class X> using Scalar = X;

template<class X> class Vector {
  public:
    typedef Vector<typename X::GenericType> GenericType;
    Vector(InitializerList<X>);
    Vector(SizeType n, std::function<X(SizeType)>);
    SizeType size() const;
    X const& operator[](SizeType i) const;
};

using RealVector = Vector<Real>;

template class Float<DoublePrecision>;
template class Bounds<Float<DoublePrecision>>;

using ScalarUnivariate = Real(Real);
using ScalarMultivariate = Real(Vector<Real>);
using VectorUnivariate = Vector<Real>(Real);
using VectorMultivariate = Vector<Real>(Vector<Real>);






template<class SIG, class X> class Differential { };

template<class X, class... AS> class Differential<Real(AS...),X> { 
    using A=Differential<Real(AS...),X>;
    
    friend A add(A, A);
    friend A mul(A, A);
    friend A add(A, X);
    friend A mul(A, X);
    friend A add(X, A);
    friend A mul(X, A);
};


template<class A> struct IsAlgebra;
template<class A, class X> struct IsAlgebraOver;

template<class X, class... AS> struct IsAlgebraOver<Differential<Real(AS...),X>,X> : True { };

template<class X> class Algebra {
    template<class A, class=EnableIf<IsAlgebraOver<A,X>>> Algebra(A const&);
    template<class A> A extract() const;
    
    friend Algebra<X> add(Algebra<X>, Algebra<X>);
    friend Algebra<X> mul(Algebra<X>, Algebra<X>);
    friend Algebra<X> add(Algebra<X>, X);
    friend Algebra<X> mul(Algebra<X>, X);
    friend Algebra<X> add(X, Algebra<X>);
    friend Algebra<X> mul(X, Algebra<X>);
};  

template<class T> using ValidatedAlgebra = Algebra<Validated<T>>;


template<template<class>class, class X> class OpenSet;

template<class X> class OpenSet<Effective,X> {
    Effective<Sierpinskian> contains(Effective<X>);
    Validated<Sierpinskian> contains(Validated<X>);
};

template<class X> class OpenSet<Validated,X> {
    Validated<Sierpinskian> contains(Validated<X>);
};

template<template<class>class P, class X> class LowerMeasurableSet;

template<class SIG> class MeasurableFunction;

template<template<class>class P, class R, class A> class MeasurableFunction<P<R(A)>> {
    LowerMeasurableSet<P,A> preimage(OpenSet<P,R>);
};

template<class SIG> class ContinuousFunction;

template<template<class>class P, class... AS> struct VectorTraits<ContinuousFunction<P<RealVector(AS...)>>> { 
    typedef ContinuousFunction<P<Real(AS...)>> ScalarType; };
//template<class... AS> struct VectorTraits<ContinuousFunction<Validated<RealVector(AS...)>>> { 
//    typedef ContinuousFunction<Validated<Real(AS...)>> ScalarType; };
//template<class... AS> struct VectorTraits<ContinuousFunction<Effective<RealVector(AS...)>>> { 
//    typedef ContinuousFunction<Effective<Real(AS...)>> ScalarType; };

template<class R, class... AS> class ContinuousFunction<Validated<R(AS...)>> 
    : public AlgebraOperations<ContinuousFunction<Validated<R(AS...)>>,Validated<R>>
{
  public:
    Validated<R> operator() (Validated<AS> const& ...);
//    ValidatedAlgebra<R> operator() (ValidatedAlgebra<AS...> const&);

    friend Validated<R> evaluate(ContinuousFunction<Validated<R(AS...)>>, Validated<AS> ...);

    template<class RR> friend auto 
    compose(ContinuousFunction<Validated<RR(R)>>, ContinuousFunction<Validated<R(AS...)>>)
        -> ContinuousFunction<Validated<RR(AS...)>>;
};
    
template<class R, class... AS> class ContinuousFunction<Effective<R(AS...)>> 
    : public AlgebraOperations<ContinuousFunction<Effective<R(AS...)>>,Effective<R>>
{
  public:
    Validated<R> operator() (Validated<AS> const& ...); 
    Effective<R> operator() (Effective<AS> const& ...); 

    friend Validated<R> evaluate(ContinuousFunction<Effective<R(AS...)>>, Validated<AS> ...);
    friend Effective<R> evaluate(ContinuousFunction<Effective<R(AS...)>>, Effective<AS> ...);
    
    template<class RR> friend auto 
        compose(ContinuousFunction<Effective<RR(R)>>, ContinuousFunction<Effective<R(AS...)>>) -> ContinuousFunction<Effective<RR(AS...)>>;

};    

using ScalarUnivariateContinuousFunction = ContinuousFunction<Real(Real)>;
template<class SIG> using EffectiveContinuousFunction = Effective<ContinuousFunction<SIG>>;
template<class SIG> using ValidatedContinuousFunction = Validated<ContinuousFunction<SIG>>;

class DomainBound { public: DomainBound(double); };
class IntervalDomain { 
  public:
    template<class X> using ElementType = Scalar<X>;
    IntervalDomain(double,double);
};
class BoxDomain { 
  public:
    template<class X> using ElementType = Vector<X>;
    BoxDomain(SizeType, Function<IntervalDomain(SizeType)>);
};

template<class... AS> struct PatchDomainTrait;
template<class... AS> using PatchDomainType = typename PatchDomainTrait<AS...>::Type;
template<> struct PatchDomainTrait<Real> { typedef IntervalDomain Type; };
template<> struct PatchDomainTrait<RealVector> { typedef BoxDomain Type; };

template<class SIG> struct SignatureTraits;
template<class SIG> using ResultType = typename SignatureTraits<SIG>::ResultType;
template<class R, class... AS> struct SignatureTraits<R(AS...)> { typedef R ResultType; };
template<template<class>class P, class R, class... AS> struct SignatureTraits<P<R(AS...)>> { typedef P<R> ResultType; };


template<class SIG> class FunctionPatch;


template<template<class>class P, class R, class... AS> class FunctionPatch<P<R(AS...)>> 
    : public AlgebraOperations<FunctionPatch<P<R(AS...)>>,P<R>>
{
  public:
    typedef P<Real> NormType;

    typedef PatchDomainType<AS...> DomainType;
    FunctionPatch(DomainType, ContinuousFunction<P<R(AS...)>>);
};

template<class SIG> FunctionPatch<SIG> restrict(ContinuousFunction<SIG>, typename FunctionPatch<SIG>::DomainType);



template<class X> struct InformationTraits;
template<class X> struct InformationTraits { typedef Validated<X> Type; };
template<class X> using EffectiveType = typename InformationTraits<Effective<X>>::Type;
template<class X> using ValidatedType = typename InformationTraits<Validated<X>>::Type;
template<class X> using ApproximateType = typename InformationTraits<Approximate<X>>::Type;
template<class X> using CanonicalType = typename InformationTraits<X>::Type;

template<class PR> struct InformationTraits<Validated<Float<PR>>> { typedef Bounds<Float<PR>> Type; };



template<class SIG, class FLT> class FunctionModel;

template<template<class>class P, class SIG, class FLT> class FunctionModel<P<SIG>,FLT>
    : public AlgebraOperations<FunctionModel<P<SIG>,FLT>,CanonicalType<P<FLT>>>
{
};
    

    

template<class X> void field_concept(X vx=X()) {
    vx=neg(vx); auto w=hlf(vx); vx=add(vx,vx); vx=sub(vx,vx); vx=mul(vx,vx); auto q=div(vx,vx); 
};

template<class A, class X> void algebra_concept(A a=A(), X vx=X()) {
    a=add(a,a); a=sub(a,a); a=mul(a,a); 
    a=add(a,vx); a=sub(a,vx); a=mul(a,vx); a=div(a,vx); 
    a=add(vx,a); a=sub(vx,a); a=mul(vx,a); 
};

template<class VA, class X> void vector_algebra_concept(VA va=VA(), X x=X()) {
    using A = decltype(declval<VA>()[declval<SizeType>()]);
    using VX = Vector<X>;
    A a=va[0]; VX vx=Vector<X>({x});
    va=add(va,va); va=sub(va,va); va=mul(a,va); va=mul(va,a); 
    va=add(va,vx); va=sub(va,vx); va=mul(va,x); va=div(va,x); 
    va=add(vx,va); va=sub(vx,va); va=mul(x,va);
};


int main() {
    Effort eff(0u);

    // Abstract logical
    Effective<Sierpinskian> l;
    Validated<Sierpinskian> vl=l.check(eff);

    // Abstract numeric
    EffectiveReal r=Rational(2,3);
    
    ValidatedReal vr(r,eff);
    vr=add(vr,vr); 
//    vr=add(vr,r);

    field_concept(r);
    field_concept(vr);
    
    DoublePrecision pr;
    Bounds<Float<DoublePrecision>> xb(0,pr);
    Approximation<Float<DoublePrecision>> xa=xb; 
    field_concept(xb);
    field_concept(xa);

    vr=xb;
#warning Would like FloatBounds vx(0,pr);

    ContinuousFunction<Effective<ScalarUnivariate>> f;
    r=f(r); r=evaluate(f,r);
    vr=f(vr); vr=evaluate(f,vr);
    f=compose(f,f);

    algebra_concept(f,r);
    
    ContinuousFunction<Validated<Real(Real,Real)>> bf;
    ContinuousFunction<Validated<Real(Real)>> uf;
   
    bf(vr,vr);
    bf=compose(uf,bf);
    
    algebra_concept(uf,vr);
    
    ContinuousFunction<Validated<RealVector(Real)>> bcf;
    vector_algebra_concept(bcf,vr);
  
}




