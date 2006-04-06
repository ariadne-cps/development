/***************************************************************************
 *            affine_vector_field.h
 *
 *  Fri Feb  4 08:57:39 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
 /*! \file affine_vector_field.h
 *  \brief vector_type fields of affine form of the form \f$\dot{x}=Ax+b\f$.
 */

#ifndef _AFFINE_VECTOR_FIELD_H
#define _AFFINE_VECTOR_FIELD_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/interval_matrix.h"

#include "../evaluation/vector_field.h"
#include "../evaluation/affine_map.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief An affine vector field in Euclidean space. */
    template <typename R>
    class AffineVectorField : public VectorField<R> 
    {
      typedef typename Geometry::Polyhedron<R> Polyhedron;
      typedef typename Geometry::Rectangle<R> Rectangle;
    
     public:
      typedef typename Geometry::Point<R> state_type;
      
      typedef LinearAlgebra::matrix<R> matrix_type;
      typedef LinearAlgebra::vector<R> vector_type;
    
      typedef LinearAlgebra::interval_vector<R> Intervalvector_type;
      typedef LinearAlgebra::interval_matrix<R> Intervalmatrix_type;
    
      virtual ~AffineVectorField();
      
      AffineVectorField(const AffineVectorField<R>& F) : _A(F.A()), _b(F.b()) { }
      AffineVectorField(const matrix_type &A, const vector_type &b) : _A(A), _b(b) { }
    
      LinearAlgebra::vector<R> apply(const Geometry::Point<R>& s) const;
      LinearAlgebra::interval_vector<R> apply(const Geometry::Rectangle<R>& r) const;
    
      LinearAlgebra::matrix<R> derivative(const Geometry::Point<R>& x) const;
      LinearAlgebra::interval_matrix<R> derivative(const Geometry::Rectangle<R>& r) const;
      
      const LinearAlgebra::matrix<R>& A() const { return this->_A; }
      const LinearAlgebra::vector<R>& b() const { return this->_b; }
      
      dimension_type dimension() const {
        return this->_b.size();
      }
      
      std::string name() const { return "AffineVectorField"; }
      
     private:
      LinearAlgebra::matrix<R> _A;
      LinearAlgebra::vector<R> _b;
    };
 
  }
}

namespace Ariadne {
  namespace LinearAlgebra {
    template <typename R>
    matrix<R> 
    exp_Ah_approx(const matrix<R>& A, 
                  const R& h, 
                  const R& e); 

    /*! \brief Compute \f$A^{-1}(e^{Ah}-I) = h\sum_{n=0}^{\infty} \frac{{(Ah)}^{n}}{(n+1)!}\f$. */
    template <typename R> 
    matrix<R> 
    exp_Ah_sub_id_div_A_approx(const matrix<R>& A, 
                               const R& h, 
                               const R& e); 

  }
}

#endif /* _AFFINE_VECTOR_FIELD_H */
