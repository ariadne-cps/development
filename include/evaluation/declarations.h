/***************************************************************************
 *            evaluation/declarations.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file evaluation/declarations.h
 *  \brief Forward declarations of classes in the Evaluation module.
 */

#ifndef ARIADNE_EVALUATION_DECLARATIONS_H
#define ARIADNE_EVALUATION_DECLARATIONS_H

namespace Ariadne { 
  

    template<class Aprx, class ES> class ApproximatorInterface;
    template<class ES> class ReducerInterface;
    template<class ES> class SubdividerInterface;
    
    template<class R> class SolverInterface;
    template<class R> class DetectorInterface;

    template<class R> class BounderInterface;
    template<class R> class FlowerInterface;

    template<class ES> class SatisfierInterface;
    template<class ES> class ApplicatorInterface;
    template<class ES, class TM> class IntegratorInterface;

    template<class Sys, class ES> class EvolverInterface;
    template<class Sys, class Aprx> class DiscretiserInterface;

    class EvolutionProfiler;
    template<class R> class EvolutionParameters;


} // namespace Ariadne

#endif /* ARIADNE_EVALUATION_DECLARATIONS_H */
