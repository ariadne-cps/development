/***************************************************************************
 *            euler_integrator.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file euler_integrator.h
 *  \brief Simple methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_EULER_INTEGRATOR_H
#define ARIADNE_EULER_INTEGRATOR_H

#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/integrator_base.h"

namespace Ariadne {
  
   

    /*! \ingroup Integrators
     *  \brief An integrator based on the Euler method.
     */
    template<class R>
    class EulerIntegrator
      : public IntegratorBase< Rectangle<R> >
    {
      typedef Interval<R> I;
     public:
      /*! \brief Constructor. */
      EulerIntegrator();


      /*! \brief Cloning operator. */
      virtual EulerIntegrator<R>* clone() const;

      /*! \brief Compute an integration time and a bounding box, given a bounding box for the intitial set, and a maximum allowable flow time. */
      virtual 
      std::pair< Rational, Box<R> >
      flow_bounds(const VectorField<R>& vector_field, 
                  const Rectangle<R>& initial_set,
                  const Rational& maximum_step_size) const; 

      /*! \brief A C0 algorithm for integrating forward a rectangle. */
      virtual Rectangle<R> 
      integration_step(const VectorField<R>& vector_field,
                       const Rectangle<R>& initial_set,
                       const Rational& maximum_step_size,
                       const Box<R>& bounding_set) const;

      /*! \brief A C0 algorithm for integrating forward a rectangle up to a certain time. */
      virtual Rectangle<R> 
      reachability_step(const VectorField<R>& vector_field,
                        const Rectangle<R>& initial_set,
                        const Rational& maximum_step_size,
                        const Box<R>& bounding_set) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
    };

    
  
} // namespace Ariadne


#endif /* ARIADNE_EULER_INTEGRATOR_H */
