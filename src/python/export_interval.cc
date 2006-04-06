/***************************************************************************
 *            python/export_interval.cc
 *
 *  21 October 2005
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

#include <boost/python.hpp>

#include "python/typedefs.h"

using namespace boost::python;
using namespace Ariadne;

/* FIXME: Do these divisions */
inline RInterval div_ii(const RInterval& r, const RInterval& i) { assert(false); }
inline RInterval div_ri(const Real& r, const RInterval& i) { assert(false); }
inline RInterval div_ir(const RInterval& i, const Real& r) { assert(false); }

void export_interval()
{
  class_< RInterval >("Interval")
    .def(init<int,int>())
    .def(init<double,double>())
    .def(init<Real,Real>())
    .def(init<double>())
    .def(init<Real>())
    .def(-self)        // __neg__
    .def(self + self)  // __add__
    .def(self - self)  // __sub__
    .def(self * self)  // __mul__
    .def("__div__",&div_ii)  // __div__
    .def(self + Real())  // __add__
    .def(self - Real())  // __sub__
    .def(self * Real())  // __mul__
    .def("__div__",&div_ir)  // __div__
    .def(Real() + self)  // __add__
    .def(Real() - self)  // __sub__
    .def(Real() * self)  // __mul__
    .def("__div__",&div_ri)  // __div__
    .def("lower", &RInterval::lower, return_value_policy<copy_const_reference>())
    .def("upper", &RInterval::upper, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

/*
  IFUN iabs(&boost::numeric::abs);
  IFUN iexp(&boost::numeric::exp);
  IFUN ilog(&boost::numeric::log);
  IFUN isin(&boost::numeric::sin);
  IFUN icos(&boost::numeric::cos);
  IFUN itan(&boost::numeric::tan);
  IFUN iasin(&boost::numeric::asin);
  IFUN iacos(&boost::numeric::acos);
  IFUN iatan(&boost::numeric::atan);
  IFUN isinh(&boost::numeric::sinh);
  IFUN icosh(&boost::numeric::cosh);
  IFUN itanh(&boost::numeric::tanh);
  IFUN iasinh(&boost::numeric::asinh);
  IFUN iacosh(&boost::numeric::acosh);
  IFUN iatanh(&boost::numeric::atanh);

  def("abs", iabs, "interval absolute value function");
  def("exp", iexp);
  def("log", ilog);
  def("sin", isin);
  def("cos", icos);
  def("tan", itan);
  def("asin", iasin);
  def("acos", iacos);
  def("atan", iatan);
  def("sinh", isinh);
  def("cosh", icosh);
  def("tanh", itanh);
  def("asinh", iasinh);
  def("acosh", iacosh);
  def("atanh", iatanh);
*/
}
