/***************************************************************************
 *            python/export_logging.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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


#include "output/logging.h"
#include "evaluation/declarations.h"

#include "python/utilities.h"
#include "python/float.h"
using namespace Ariadne;
using namespace Ariadne::Output;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;


void export_logging()
{
  def("set_linear_algebra_verbosity",&set_linear_algebra_verbosity);
  def("set_combinatoric_verbosity",&set_combinatoric_verbosity);
  def("set_geometry_verbosity",&set_geometry_verbosity);
  def("set_evaluation_verbosity",&set_evaluation_verbosity);
  def("set_applicator_verbosity",&set_applicator_verbosity);
  def("set_integrator_verbosity",&set_integrator_verbosity);
  def("set_hybrid_evolver_verbosity",&set_hybrid_evolver_verbosity);
  def("set_input_verbosity",&set_input_verbosity);

  def("redirect_log",&redirect_log);
}

