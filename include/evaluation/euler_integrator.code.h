/***************************************************************************
 *            euler_integrator.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
//#define DEBUG

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "euler_integrator.h"

#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"

#include "../system/vector_field.h"

#include "../output/logging.h"

namespace Ariadne {
    

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R>
Evaluation::EulerIntegrator<R>::EulerIntegrator()
{
}



template<class R>
Evaluation::EulerIntegrator<R>*
Evaluation::EulerIntegrator<R>::clone() const
{
  return new EulerIntegrator<R>();
}



template<class R>
Geometry::Point< Numeric::Interval<R> >
Evaluation::EulerIntegrator<R>::flow_step(const System::VectorFieldInterface<R>& vector_field, 
                                          const Geometry::Point<I>& initial_point, 
                                          const Numeric::Interval<R>& step_size, 
                                          const Geometry::Rectangle<R>& bounding_set) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"EulerIntegrator::flow_step(VectorFieldInterface,Point,Interval,Rectangle) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_point,"EulerIntegrator::flow_step(VectorFieldInterface,Point,Interval,Rectangle) const");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,bounding_set,"EulerIntegrator::flow_step(VectorFieldInterface,Point,Interval,Rectangle) const");
  
  return initial_point + ( I(step_size) * vector_field(Point<I>(bounding_set)) );
}





template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                         const Geometry::Rectangle<R>& initial_set, 
                                                         const Numeric::Interval<R>& step_size, 
                                                         const Geometry::Rectangle<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,Interval,Rectangle) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,Interval,Rectangle) const");

  return initial_set + I(step_size) * vector_field(bounding_set);
}



template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                          const Geometry::Rectangle<R>& initial_set, 
                                                          const Numeric::Interval<R>& step_size, 
                                                          const Geometry::Rectangle<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,Numeric::Interval<R>) const\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set(),"EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,Numeric::Interval<R>) const");
  
  return initial_set + I(0,step_size.upper()) * vector_field(bounding_set);
}






template<class R>
Geometry::Zonotope<typename Evaluation::EulerIntegrator<R>::I,R> 
Evaluation::EulerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                 const Geometry::Zonotope<I,R>& initial_set, 
                                                 const Numeric::Interval<R>& step_size, 
                                                 const Geometry::Rectangle<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::integration_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::integration_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const");
  return Geometry::Zonotope<I,R>(this->integration_step(vector_field,initial_set.bounding_box(),step_size,bounding_set));
}




template<class R>
Geometry::Zonotope<typename Evaluation::EulerIntegrator<R>::I,R> 
Evaluation::EulerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,R>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::reachability_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::reachability_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const");
  return Geometry::Zonotope<I,R>(this->reachability_step(vector_field,initial_set.bounding_box(),step_size,bounding_set));
}



template<class R>
Geometry::Zonotope<typename Evaluation::EulerIntegrator<R>::I> 
Evaluation::EulerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                 const Geometry::Zonotope<I,I>& initial_set, 
                                                 const Numeric::Interval<R>& step_size, 
                                                 const Geometry::Rectangle<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::integration_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::integration_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const");
  return Geometry::Zonotope<I,I>(this->integration_step(vector_field,initial_set.bounding_box(),step_size,bounding_set));
}




template<class R>
Geometry::Zonotope<typename Evaluation::EulerIntegrator<R>::I> 
Evaluation::EulerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,I>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::reachability_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::reachability_step(VectorFieldInterface,Zonotope,Interval,Rectangle) const");
  return Geometry::Zonotope<I,I>(this->reachability_step(vector_field,initial_set.bounding_box(),step_size,bounding_set));
}




}

