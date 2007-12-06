/***************************************************************************
 *            python/export_set_based_hybrid_evolver.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
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

#include "python/float.h"

#include "geometry/zonotope.h"
#include "geometry/hybrid_set.h"
#include "system/set_based_hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/set_based_hybrid_evolver.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_set_based_hybrid_evolver() 
{

  typedef Numeric::Interval<R> I;
  typedef typename Evaluation::VectorFieldEvolver<R>::basic_set_type BS;

  class_< SetBasedHybridEvolver<R> > evolver_class("SetBasedHybridEvolver",init< EvolutionParameters<R> >());
  evolver_class.def(init<const EvolutionParameters<R>&,const ApplicatorInterface<BS>&,const IntegratorInterface<BS>&>());
  evolver_class.def(init<const SetBasedHybridEvolver<R>&>());
  evolver_class.def("discrete_step",
                    (HybridSet<R>(SetBasedHybridEvolver<R>::*)(const SetBasedHybridAutomaton<R>&,const HybridSet<R>&))&SetBasedHybridEvolver<R>::discrete_step);
  evolver_class.def("continuous_chainreach",
                    (HybridSet<R>(SetBasedHybridEvolver<R>::*)(const SetBasedHybridAutomaton<R>&,const HybridSet<R>&))&SetBasedHybridEvolver<R>::continuous_chainreach);
  evolver_class.def("chainreach",
                    (HybridSet<R>(SetBasedHybridEvolver<R>::*)(const SetBasedHybridAutomaton<R>&,const HybridSet<R>&))&SetBasedHybridEvolver<R>::chainreach);


  evolver_class.def("continuous_chainreach",
                    (HybridSet<R>(SetBasedHybridEvolver<R>::*)(const SetBasedHybridAutomaton<R>&,const HybridSet<R>&,const HybridSet<R>&))&SetBasedHybridEvolver<R>::continuous_chainreach);
  evolver_class.def("chainreach",
                    (HybridSet<R>(SetBasedHybridEvolver<R>::*)(const SetBasedHybridAutomaton<R>&,const HybridSet<R>&,const HybridSet<R>&))&SetBasedHybridEvolver<R>::chainreach);
  
}


template void export_set_based_hybrid_evolver<FloatPy>();
