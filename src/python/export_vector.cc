/***************************************************************************
 *            python/export_vector.cc
 *
 *  17 November 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can rediself_ns::stribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "base/numerical_type.h"
#include "base/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/interval_vector.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

inline Real rvector_getitem(const RVector& v, uint i) {
  return v(i);
}

inline void rvector_setitem(RVector& v, uint i, Real x) {
  v(i)=x;
}

inline void rvector_setitem_from_double(RVector& v, uint i, double x) {
  v(i)=Ariadne::convert_to<Real>(x);
}

inline RVector rvector_add_rvector(const RVector& u, RVector& v) {
  return RVector(u+v);
}

inline RVector rvector_sub_rvector(const RVector& u, RVector& v) {
  return RVector(u-v);
}

inline RInterval ivector_getitem(const RIntervalVector& v, uint i) {
  return v(i);
}

inline void ivector_setitem(RIntervalVector& v, uint i, RInterval x) {
  v(i)=x;
}

inline Field qvector_getitem(const FVector& v, uint i) {
  return v(i);
}

inline void qvector_setitem(FVector& v, uint i, Field x) {
  v(i)=x;
}

inline void qvector_setitem_from_real(FVector& v, uint i, Real x) {
  v(i)=Ariadne::convert_to<Rational>(x);
}

inline void qvector_setitem_from_double(FVector& v, uint i, double x) {
  v(i)=Ariadne::convert_to<Rational>(x);
}


void export_vector() {
  class_<RVector>("Vector",init<int>())
    .def(init<RVector>())
    .def("__len__", &RVector::size)
    .def("__getitem__",&rvector_getitem)
    .def("__setitem__",&rvector_setitem)
    .def("__setitem__",&rvector_setitem_from_double)
    .def("__add__",&rvector_add_rvector)
    .def("__sub__",&rvector_sub_rvector)
    .def(self_ns::str(self))    // __self_ns::str__
  ;

  class_<FVector>("KVector",init<int>())
    .def(init<FVector>())
    .def("__len__", &FVector::size)
    .def("__getitem__",&qvector_getitem)
    .def("__setitem__",&qvector_setitem)
    .def("__setitem__",&qvector_setitem_from_real)
    .def("__setitem__",&qvector_setitem_from_double)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
}

void export_interval_vector() {
  class_<RIntervalVector>("IntervalVector",init<int>())
    .def(init<RIntervalVector>())
    .def("__len__", &RIntervalVector::size)
    .def("__getitem__",&ivector_getitem)
    .def("__setitem__",&ivector_setitem)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
