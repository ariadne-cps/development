/***************************************************************************
 *            taylor_series.inline.h
 *
 *  Copyright 2007  A Pieter Collins
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
 
#include "linear_algebra/vector.h"
#include "function/function_series.h"

namespace Ariadne {

template<class X> inline
Function::TaylorSeries<X>::TaylorSeries() 
  : _data(1u)
{
}

template<class X> inline
Function::TaylorSeries<X>::TaylorSeries(smoothness_type d) 
  : _data(d+1u)
{
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>::TaylorSeries(smoothness_type d, const XX* ptr) 
  : _data(ptr,ptr+d+1)
{
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>::TaylorSeries(const TaylorSeries<XX>& ts) 
  : _data(ts.data())
{
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>&
Function::TaylorSeries<X>::operator=(const TaylorSeries<XX>& ts) 
{
  this->_data=ts.data();
  return *this;
}

template<class X> template<class XX> inline
Function::TaylorSeries<X>&
Function::TaylorSeries<X>::operator=(const XX& c) 
{
  this->_data[0]=c;
  return *this;
}

template<class X> inline
smoothness_type 
Function::TaylorSeries<X>::degree() const 
{
  return this->_data.size()-1;
}

template<class X> inline 
const X& 
Function::TaylorSeries<X>::value() const
{
  return this->_data[0];
}

template<class X> inline 
X& 
Function::TaylorSeries<X>::value() 
{
  return this->_data[0];
}

template<class X> inline
const array<X>& 
Function::TaylorSeries<X>::data() const
{
  return this->_data;
}

template<class X> inline
array<X>& 
Function::TaylorSeries<X>::data() 
{
  return this->_data;
}

template<class X> inline 
const X& 
Function::TaylorSeries<X>::operator[](const smoothness_type& j) const
{
  return this->_data[j];
}

template<class X> inline
X& 
Function::TaylorSeries<X>::operator[](const smoothness_type& j)
{
  return this->_data[j];
}




template<class X> inline 
Function::TaylorSeries<X> 
Function::min(const TaylorSeries<X>& x1, const TaylorSeries<X>& x2) 
{
  if(x1[0]==x2[0]) {
    ARIADNE_THROW(std::runtime_error,"min(TaylorSeries x1, TaylorSeries x2)","x1[0]==x2[0]");
  }
  return x1[0]<x2[0] ? x1 : x2;
}

template<class X> inline 
Function::TaylorSeries<X> 
Function::max(const TaylorSeries<X>& x1,const TaylorSeries<X>& x2) 
{
  if(x1[0]==x2[0]) { 
    ARIADNE_THROW(std::runtime_error,"max(TaylorSeries x1, TaylorSeries x2)","x1[0]==x2[0]"); 
  }
  return x1[0]>x2[0] ? x1 : x2;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::pos(const TaylorSeries<X>& x)
{
  return x;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::neg(const TaylorSeries<X>& x)
{
  TaylorSeries<X> result(x.degree());
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = -x[n];
  }
  return result;
}

template<class X> inline 
Function::TaylorSeries<X> 
Function::abs(const TaylorSeries<X>& x) 
{
  if(x[0]==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(TaylorSeries x)","x[0]==0"); 
  }
  return x[0]>0 ? pos(x) : neg(x); 
}

template<class X> inline
Function::TaylorSeries<X> 
Function::rec(const TaylorSeries<X>& x)
{
  return compose(FunctionSeries<X>::rec(x.degree(),x.value()),x);
}



template<class X> inline
Function::TaylorSeries<X> 
Function::add(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]+y[n];
  }
  return result;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::sub(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n] = x[n]-y[n];
  }
  return result;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::mul(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n]=x[0]*y[n];
    for(size_type i=1; i<=n; ++i) {
      result[n] += Numeric::bin<int>(n,i)*x[i]*y[n-i];
    }
  }
  return result;
}

template<class X> inline
Function::TaylorSeries<X> 
Function::div(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return x*rec(y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::pow(const TaylorSeries<X>& x, const uint& k)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),k),x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::pow(const TaylorSeries<X>& x, const int& k)
{
  return compose(FunctionSeries<X>::pow(x.degree(),x.value(),uint(k)),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::sqrt(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::exp(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::exp(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::log(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::log(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::sin(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::sin(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::cos(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::cos(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::tan(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::tan(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::asin(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::asin(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::acos(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::acos(x.degree(),x.value()),x);
}

template<class X>  
Function::TaylorSeries<X> 
Function::atan(const TaylorSeries<X>& x) 
{
  return compose(FunctionSeries<X>::atan(x.degree(),x.value()),x);
}


template<class X> inline
Function::TaylorSeries<X> 
Function::operator-(const TaylorSeries<X>& x)
{
  return neg(x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator+(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return add(x,y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator-(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return sub(x,y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator*(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return mul(x,y);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator/(const TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  return div(x,y);
}


template<class X> inline
Function::TaylorSeries<X>& 
Function::operator+=(TaylorSeries<X>& x, const TaylorSeries<X>& y)
{
  for(size_type n=0; n<=std::min(x.degree(),y.degree()); ++n) {
    x[n] += y[n];
  }
  return x;
}



template<class X> inline
Function::TaylorSeries<X> 
Function::operator+(const TaylorSeries<X>& x, const X& c)
{
  return add(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator+(const X& c, const TaylorSeries<X>& x)
{
  return add(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator-(const TaylorSeries<X>& x, const X& c)
{
  return sub(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator-(const X& c, const TaylorSeries<X>& x)
{
  return sub(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator*(const TaylorSeries<X>& x, const X& c)
{
  return mul(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator*(const X& c, const TaylorSeries<X>& x)
{
  return mul(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator/(const TaylorSeries<X>& x, const X& c)
{
  return div(x,TaylorSeries<X>::constant(x.degree(),c));
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator/(const X& c, const TaylorSeries<X>& x)
{
  return div(TaylorSeries<X>::constant(x.degree(),c),x);
}

template<class X> inline
Function::TaylorSeries<X> 
Function::operator/(const double& c, const TaylorSeries<X>& x)
{
  return X(c)/x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator+=(TaylorSeries<X>& x, const X& c)
{
  x[0]+=c;
  return x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator-=(TaylorSeries<X>& x, const X& c)
{
  x[0]-=c;
  return x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator*=(TaylorSeries<X>& x, const X& c)
{
  reinterpret_cast< LinearAlgebra::Vector<X>& >(x.data())*=c;
  return x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator/=(TaylorSeries<X>& x, const X& c)
{
  reinterpret_cast< LinearAlgebra::Vector<X>& >(x.data())/=c;
  return x;
}


template<class X>  
Function::TaylorSeries<X>&
Function::operator+=(TaylorSeries<X>& x, const double& c)
{
  X& v=x[0];
  v+=c;
  return x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator-=(TaylorSeries<X>& x, const double& c)
{
  x[0]-=c;
  return x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator*=(TaylorSeries<X>& x, const double& c)
{
  reinterpret_cast< LinearAlgebra::Vector<X>& >(x.data())*=c;
  return x;
}

template<class X>  
Function::TaylorSeries<X>&
Function::operator/=(TaylorSeries<X>& x, const double& c)
{
  reinterpret_cast< LinearAlgebra::Vector<X>& >(x.data())/=c;
  return x;
}


template<class X, class R>  
Function::TaylorSeries<X>&
Function::operator*=(TaylorSeries<X>& x, const R& c)
{
  reinterpret_cast< LinearAlgebra::Vector<X>& >(x.data())*=c;
  return x;
}







} // namespace Ariadne
