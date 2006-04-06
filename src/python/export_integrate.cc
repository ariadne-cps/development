/***************************************************************************
 *            python/export_integrate.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"

#include "evaluation/integrate.h"


#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

RRectangle rectangle_integration_step(const RVectorFieldBase& vf, const RRectangle& r, Real& h) {
  return integration_step(vf,r,h);
}

void export_integrate() {
  typedef RRectangle (*IntStepRectFunc) (const RVectorFieldBase&, const RRectangle&, Real&);
  typedef RParallelotope (*IntStepPltpFunc) (const RVectorFieldBase&, const RParallelotope&, Real&);
  typedef RZonotope (*IntStepZntpFunc) (const RVectorFieldBase&, const RZonotope&, Real&);
  typedef RRectangle (*IntRectFunc) (const RVectorFieldBase&, const RRectangle&, const Real&, const Real&);
  typedef RParallelotope (*IntPltpFunc) (const RVectorFieldBase&, const RParallelotope&, const Real&, const Real&);
  typedef RParallelotopeListSet (*IntLSPltpFunc) (const RVectorFieldBase&, const RParallelotopeListSet&, const Real&, const Real&);
  typedef RGridMaskSet (*IntGMSFunc) (const RVectorFieldBase&, const RGridMaskSet&, const Real&, const Real&);
  
  
  def("integration_step", IntStepRectFunc(&integration_step), "integrate a vector field over a set for a time up to time h");
  def("integration_step", IntStepPltpFunc(&integration_step));
  def("reach_step", IntStepRectFunc(&reach_step));
  def("reach_step", IntStepPltpFunc(&reach_step));
  def("reach_step", IntStepZntpFunc(&reach_step));
  def("integrate", IntRectFunc(&integrate));
  def("integrate", IntPltpFunc(&integrate));
  def("integrate", IntGMSFunc(&integrate));
  def("reach", IntGMSFunc(&reach));
//  def("chainreach", &chainreach_of_rectangle_list_set, "chain reach of a set" );
}
