/***************************************************************************
 *            numerical_traits.h
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
 
/*! \file numerical_traits.h
 *  \brief Traits classes to define properties of numerical types.
 */

#ifndef _ARIADNE_NUMERICAL_TRAITS_H
#define _ARIADNE_NUMERICAL_TRAITS_H

#include <string>

#include "../declarations.h"

namespace Ariadne {
  namespace Numeric {
    class Float64;
    class MPFloat;
    class Rational;
    template<class R> class Interval;
      
    /* numerical traits */
    /*! \brief Tags a class representing a ring. */
    class ring_tag { };
    /*! \brief Tags a class representing a field. */
    class field_tag { };
      
    /*! \brief Typedef's describing a numerical type. */
    template<class T1, class T2> class traits { };
  //    template<class T1, class T2=T1> class traits { };
    
    template<> struct traits<int> {
      typedef int arithmetic_type;
    };
  
    template<> struct traits<double> {
      typedef double arithmetic_type;
      typedef Interval<double> interval_type;
    };
  
    template<> struct traits< Float64 > { 
      typedef Interval<Float64> arithmetic_type; 
      typedef Interval<Float64> interval_type;
    };
    
    template<> struct traits< MPFloat > { 
      typedef Interval<MPFloat> arithmetic_type; 
      typedef Interval<MPFloat> interval_type; 
    };
    
    template<> struct traits< Rational > { 
      typedef Rational arithmetic_type; 
      typedef Interval<Rational> interval_type; 
    };

    template<class R> struct traits< Interval<R> > { 
      typedef Interval<R> arithmetic_type; 
      typedef Interval<R> interval_type; 
    };

    
    template<> struct traits< Float64, double > { 
      typedef Interval<Float64> arithmetic_type; 
    };
    
    template<> struct traits< double, Float64 > { 
      typedef Interval<Float64> arithmetic_type; 
    };
    
    template<> struct traits< Float64, Rational > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Rational, Float64 > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< MPFloat, Rational > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Rational, MPFloat > { 
      typedef Rational arithmetic_type; 
    };
    

    template<class R> struct traits< R,Interval<R> > { 
      typedef Interval<R> arithmetic_type; 
    };
    
    template<class R> struct traits< Interval<R>, R > { 
      typedef Interval<R> arithmetic_type; 
    };    
    



    //! \name Numerical type description.
    //@{
    /*! \brief The name of class T. */
    template<class T> inline std::string name();

    //! \name Standard conversion operations. (Deprecated) 
    /*! \brief Approximate \a x by an element of Res. */
    template<class Res, class Arg> inline Res convert_to(const Arg& x) { return Res(x); }
    
    /*! \brief Approximate \a x by an element of Res with accuracy \a e. */
    template<class Res, class Arg, class Err> Res approximate(const Arg& x, const Err& e);
    //@}
  }   
}
  

#endif /* _ARIADNE_NUMERICAL_TRAITS_H */
