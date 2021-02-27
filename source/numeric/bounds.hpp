/***************************************************************************
 *            numeric/bounds.hpp
 *
 *  Copyright  2008-21  Pieter Collins
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

/*! \file numeric/bounds.hpp
 *  \brief Bounds on numbers with algebraic endpoints.
 */

#ifndef ARIADNE_BOUNDS_HPP
#define ARIADNE_BOUNDS_HPP

namespace Ariadne {


template<class F> class Bounds;
template<class F, class FE> class Ball;

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Bounds<Y>;
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Ball<Y>;

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Bounds<Y> {
    Y _l, _u;
  public:
    Bounds<Y>(Y w) : _l(w), _u(w) { }
    Bounds<Y>(Y l, Y u) : _l(l), _u(u) { }
    template<class X> requires Constructible<Y,X> Bounds<Y>(Bounds<X> const& x)
        : Bounds<Y>(Y(x.lower_raw()),Y(x.upper_raw())) { }
    operator ValidatedNumber() const;
    Bounds<Y> pm(Y e) { return Bounds<Y>(_l-e,_u+e); }
    Y lower() const { return _l; }
    Y upper() const { return _u; }
    Y lower_raw() const { return _l; }
    Y upper_raw() const { return _u; }
    Bounds<FloatDP> get(DP) const;
    Bounds<FloatMP> get(MP) const;
    friend Bounds<Y> operator+(Bounds<Y> const& w) { return Bounds<Y>(+w._l,w._u); }
    friend Bounds<Y> operator-(Bounds<Y> const& w) { return Bounds<Y>(-w._u,-w._l); }
    friend Bounds<Y> operator+(Bounds<Y> const& w1, Bounds<Y> const& w2) { return Bounds<Y>(w1._l+w2._l,w1._u+w2._u); }
    friend Bounds<Y> operator-(Bounds<Y> const& w1, Bounds<Y> const& w2) { return Bounds<Y>(w1._l-w2._u,w1._u-w2._l); }
    friend Bounds<Y> operator*(Bounds<Y> const& w1, Bounds<Y> const& w2) { return Bounds<Y>::_mul(w1,w2); }
    friend Bounds<Y> operator/(Bounds<Y> const& w1, Bounds<Y> const& w2) { return Bounds<Y>::_div(w1,w2); }
    friend Bounds<Y> sqr(Bounds<Y> const& w) { return Bounds<Y>::_sqr(w); }
    friend Bounds<Y> hlf(Bounds<Y> const& w) { return Bounds<Y>(hlf(w._l),hlf(w._u)); }
    friend Bounds<Y> rec(Bounds<Y> const& w) { return Bounds<Y>::_rec(w); }
    friend Bounds<Y> abs(Bounds<Y> const& w) { return Bounds<Y>(max(min(w._l,-w._u),0),max(-w._l,w._u)); }
    friend Bounds<Y> max(Bounds<Y> const& w1, Bounds<Y> const& w2) { return Bounds<Y>(max(w1._l,w2._l),max(w1._u,w2._u)); }
    friend Bounds<Y> min(Bounds<Y> const& w1, Bounds<Y> const& w2) { return Bounds<Y>(min(w1._l,w2._l),min(w1._u,w2._u)); }
    friend ValidatedKleenean sgn(Bounds<Y> const& w) {
        if (w._l>0) { return true; } else if (w._u<0) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator==(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        if (w1._l>=w2._u && w1._u<=w2._l) { return true; } else if (w1._u< w2._l || w1._l> w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator!=(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        if (w1._u< w2._l || w1._l> w2._u) { return true; } else if (w1._l>=w2._u && w1._u<=w2._l) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator<=(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        if (w1._u<=w2._l) { return true; } else if (w1._l> w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator>=(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        if (w1._l>=w2._u) { return true; } else if (w1._u< w2._l) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator< (Bounds<Y> const& w1, Bounds<Y> const& w2) {
        if (w1._u< w2._l) { return true; } else if (w1._l>=w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator> (Bounds<Y> const& w1, Bounds<Y> const& w2) {
        if (w1._l> w2._u) { return true; } else if (w1._u<=w2._l) { return false; } else { return indeterminate; } }
    friend Bool consistent(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        return w1._l<=w2._u && w1._u >= w2._l; }
    friend Boolean refines(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        return w1._l>=w2._l and w1._u<=w2._u; }
    friend Bounds<Y> refinement(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        return Bounds<Y>(max(w1._l,w2._l),min(w1._u,w2._u)); }
    friend Bounds<Y> coarsening(Bounds<Y> const& w1, Bounds<Y> const& w2) {
        return Bounds<Y>(min(w1._l,w2._l),max(w1._u,w2._u)); }
    friend OutputStream& operator<<(OutputStream& os, Bounds<Y> y) {
        return os << "[" << y._l << ":" << y._u << "]"; }
  private:
    static Bounds<Y> _mul(Bounds<Y> const& y1, Bounds<Y> const& y2) {
        using B=Bounds<Y>; const Y& y1l=y1._l; const Y& y1u=y1._u; const Y& y2l=y2._l; const Y& y2u=y2._u;
        if(y1l>=0) {
            if(y2l>=0) { return B(y1l*y2l,y1u*y2u); } else if(y2u<=0) { return B(y1u*y2l,y1l*y2u); } else { return B(y1u*y2l,y1u*y2u); } }
        else if(y1u<=0) {
            if(y2l>=0) { return B(y1l*y2u,y1u*y2l); } else if(y2u<=0) { return B(y1u*y2u,y1l*y2l); } else { return B(y1l*y2u,y1l*y2l); } }
        else {
            if(y2l>=0) { return B(y1l*y2u,y1u*y2u); } else if(y2u<=0) { return B(y1u*y2l,y1l*y2l); }
            else { return B(min(y1u*y2l,y1l*y2u),max(y1l*y2l,y1u*y2u)); } } }
    static Bounds<Y> _sqr(Bounds<Y> const& y) {
        if(y._l>0) { return Bounds<Y>(sqr(y._l),sqr(y._u)); }
        else if(y._u<0) { return Bounds<Y>(sqr(y._u),sqr(y._l)); }
        else { return Bounds<Y>(nul(y._l),max(sqr(y._l),sqr(y._u))); } }

    static Bounds<Y> _div(Bounds<Y> const& y1, Bounds<Y> const& y2);
    static Bounds<Y> _rec(Bounds<Y> const& y);
};

/*
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y>
Bounds<Y> Bounds<Y>::_mul(Bounds<Y> const& y1, Bounds<Y> const& y2) {
    const Y& y1l=y1.lower_raw(); const Y& y1u=y1.upper_raw();
    const Y& y2l=y2.lower_raw(); const Y& y2u=y2.upper_raw();
    if(y1l>=0) {
        if(y2l>=0) {
            return Bounds<Y>(mul(y1l,y2l),mul(y1u,y2u));
        } else if(y2u<=0) {
            return Bounds<Y>(mul(y1u,y2l),mul(y1l,y2u));
        } else {
            return Bounds<Y>(mul(y1u,y2l),mul(y1u,y2u));
        }
    }
    else if(y1u<=0) {
        if(y2l>=0) {
            return Bounds<Y>(mul(y1l,y2u),mul(y1u,y2l));
        } else if(y2u<=0) {
            return Bounds<Y>(mul(y1u,y2u),mul(y1l,y2l));
        } else {
            return Bounds<Y>(mul(y1l,y2u),mul(y1l,y2l));
        }
    } else {
        if(y2l>=0) {
            return Bounds<Y>(mul(y1l,y2u),mul(y1u,y2u));
        } else if(y2u<=0) {
            return Bounds<Y>(mul(y1u,y2l),mul(y1l,y2l));
        } else {
            return Bounds<Y>(min(mul(y1u,y2l),mul(y1l,y2u)),max(mul(y1l,y2l),mul(y1u,y2u)));
        }
    }
}

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> auto
Bounds<Y>::_sqr(Bounds<Y> const& y) -> Bounds<Y> {
    if(y._l>0) { return Bounds<Y>(sqr(y._l),sqr(y._u)); }
    else if(y._u<0) { return Bounds<Y>(sqr(y._u),sqr(y._l)); }
    else { return Bounds<Y>(nul(y._l),max(sqr(y._l),sqr(y._u))); }
}
*/


template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Ball<Y> {
    Y _v, _e;
  public:
    Ball<Y,Y>(Y w) : _v(w), _e(0) { }
    Ball<Y,Y>(Y v, Y e) : _v(v), _e(e) { }
    template<class X, class XE> requires Constructible<Y,X> and Constructible<Y,XE>
    Ball<Y,Y>(Ball<X,XE> const& x)
        : Ball<Y,Y>(Y(x.value_raw()),Y(x.error_raw())) { }
    explicit Ball<Y,Y>(Bounds<Y> const& w) : _v(hlf(w.lower()+w.upper())), _e(hlf(w.upper()-w.lower())) { }
    explicit operator Bounds<Y>() const { return Bounds<Y>(_v-_e,_v+_e); }
    operator ValidatedNumber() const;
    Y value() const { return _v; }
    Y error() const { return _e; }
    Y value_raw() const { return _v; }
    Y error_raw() const { return _e; }
    friend Ball<Y> operator+(Ball<Y> const& w1, Ball<Y> const& w2) { return Ball<Y>(w1._v+w2._v,w1._e+w2._e); }
    friend Ball<Y> operator-(Ball<Y> const& w1, Ball<Y> const& w2) { return Ball<Y>(w1._v-w2._v,w1._e+w2._e); }
    friend Ball<Y> operator*(Ball<Y> const& w1, Ball<Y> const& w2) { return Ball<Y>(w1._v*w2._v,abs(w1._v)*w2._e+w1._e*abs(w2._v)+w1._e*w2._e); }
    friend Ball<Y> abs(Ball<Y> const& w) {
        if (abs(w._v)>=w._e) { return Ball<Y>(abs(w._v),w._e); } else { Y av=hlf(max(w._e-w._v,w._v+w._e)); return Ball<Y>(av,av); } }
    friend ValidatedKleenean operator<(Ball<Y> const& w1, Ball<Y> const& w2) {
        if (w1._v+w1._e<w2._v-w2._e) { return true; } else if (w1._v-w1._e >= w2._v+w2._e) { return false; } else { return indeterminate; } }
    friend Boolean refines(Ball<Y> const& w1, Ball<Y> const& w2) { return abs(w1._v-w2._v)+w1._e <= w2._e; }
    friend Ball<Y> refinement(Ball<Y> const& w1, Ball<Y> const& w2) { return Ball<Y>(refinement(Bounds<Y>(w1),Bounds<Y>(w2))); }
    friend OutputStream& operator<<(OutputStream& os, Ball<Y> y) { return os << "[" << y._v << ":" << y._e << "]"; }
};

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class LowerBound<Y> {
    Y _l;
  public:
    LowerBound<Y>(Y l) : _l(l) { }
    LowerBound<Y>(Bounds<Y> lu) : _l(lu.lower_raw()) { }
    template<class X> requires Constructible<Y,X>
        LowerBound<Y>(LowerBound<X> const& x) : LowerBound<Y>(Y(x.raw())) { }
    Y raw() const { return _l; }
    LowerBound<FloatDP> get(DoublePrecision pr) const;
    LowerBound<FloatMP> get(MultiplePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, LowerBound<Y> const& y) { return os << y._l << ":"; }
};

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class UpperBound<Y> {
    Y _u;
  public:
    explicit UpperBound<Y>(Y u) : _u(u) { }
    UpperBound<Y>(Bounds<Y> lu) : _u(lu.upper_raw()) { }
    template<class X> requires Constructible<Y,X>
        UpperBound<Y>(UpperBound<X> const& x) : UpperBound<Y>(Y(x.raw())) { }
    Y raw() const { return _u; }
    UpperBound<FloatDP> get(DoublePrecision pr) const;
    UpperBound<FloatMP> get(MultiplePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, UpperBound<Y> const& y) { return os << ":" << y._u; }
};

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Approximation<Y> {
    Y _a;
  public:
    Approximation<Y>(Y a) : _a(a) { }
    template<class X> requires Constructible<Y,X>
        Approximation<Y>(Approximation<X> const& x) : Approximation<Y>(Y(x.raw())) { }
    Y raw() const { return _a; }
    Approximation<FloatDP> get(DoublePrecision pr) const;
    Approximation<FloatMP> get(MultiplePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, Approximation<Y> const& y) { return os << "~" << y._a; }
};

template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Positive<Bounds<Y>> : public Bounds<Y> {
    public: Positive(Bounds<Y> w) : Bounds<Y>(w) { };
};
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Positive<LowerBound<Y>> : public LowerBound<Y> {
    public: Positive(LowerBound<Y> w) : LowerBound<Y>(w) { };
};
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Positive<UpperBound<Y>> : public UpperBound<Y> {
    public: Positive(UpperBound<Y> w) : UpperBound<Y>(w) { } ;
};
template<class Y> requires (not RoundedRing<Y>) and LatticeRing<Y> class Positive<Approximation<Y>> : public Approximation<Y> {
    public: Positive(Approximation<Y> w) : Approximation<Y>(w) { }
};




} // namespace Ariadne

#endif
