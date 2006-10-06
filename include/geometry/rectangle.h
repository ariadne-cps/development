/***************************************************************************
 *            rectangle.h
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file rectangle.h
 *  \brief Rectangles and cuboids.
 */

#ifndef _ARIADNE_RECTANGLE_H
#define _ARIADNE_RECTANGLE_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/array.h"

#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

#include "../geometry/ppl_polyhedron.h"
#include "../geometry/point.h"
#include "../geometry/rectangle_expression.h"

namespace Ariadne {
  namespace Geometry {

    template<> 
    inline bool is_a<Rectangle,Rectangle>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Parallelotope>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Zonotope>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Polyhedron>() { return true; }

    /*! \ingroup BasicSet
     *  \brief A cuboid of arbitrary dimension.
     */
    template <typename R>
    class Rectangle 
      : public RectangleExpression< Rectangle<R> >
    {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> state_type;
    
     public:
      //@{
      //! \name Constructors
      /*! \brief Construct an empty rectangle with dimension \a d. */
      explicit Rectangle(size_type d=0)
        : _bounds(2*d)
      { 
        if(d!=0) { this->_bounds[0]=1; this->_bounds[1]=0; }
      }
      
      /*! \brief Construct a rectangle from a range of interval values. */
      template<class ForwardIterator>
      Rectangle(ForwardIterator b, ForwardIterator e)
        : _bounds(2*std::distance(b,e))
      {
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,b->lower());
          this->set_upper_bound(i,b->upper());
          ++b;
        }
      }
      
      /*! \brief Construct from an array of intervals. */
      explicit Rectangle(const array< Interval<R> >& a)
        : _bounds(2*a.size())
      {
        for(dimension_type i=0; i!=a.size(); ++i) {
          this->set_lower_bound(i,a[i].lower());
          this->set_upper_bound(i,a[i].upper());
        }
      }
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Rectangle(const std::vector< Interval<R> >& v)
        : _bounds(2*v.size())
      {
        for(dimension_type i=0; i!=v.size(); ++i) {
          this->set_lower_bound(i,v[i].lower());
          this->set_upper_bound(i,v[i].upper());
        }
      }

      /*! \brief Construct a degenerate rectangle from a single point. */
      explicit Rectangle(const Point<R>& pt)
        : _bounds(2*pt.dimension())
      {
        for(dimension_type i=0; i!=pt.dimension(); ++i) {
          this->set_lower_bound(i,pt[i]);
          this->set_upper_bound(i,pt[i]);
        }
      }
      
      /*! \brief Construct from two corners. */
      explicit Rectangle(const Point<R>& pt1, const Point<R>& pt2) 
        : _bounds(2*pt1.dimension())
      {
        if (pt1.dimension()!=pt2.dimension()) {
          throw std::domain_error("The parameters have different space dimensions");
        }
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,Numeric::min(pt1[i],pt2[i]));
          this->set_upper_bound(i,Numeric::max(pt1[i],pt2[i]));
        }
      }
      
      /*! \brief Construct from a string literal. */
      explicit Rectangle(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Rectangle(const LinearAlgebra::Vector< Interval<R> >& iv)
        : _bounds(2*iv.size())
      {
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,iv[i].lower());
          this->set_upper_bound(i,iv[i].upper());
        }
      }

      /*! \brief Convert from a rectangle expression. */
      template<class E>
      Rectangle(const RectangleExpression<E>& original)
        : _bounds(2*original().dimension())
      {         
        const E& expression=original();
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->_bounds[2*i]=expression.lower_bound(i);
          this->_bounds[2*i+1]=expression.upper_bound(i);
        }
      }
    
    
      /*! \brief Copy constructor. */
      Rectangle(const Rectangle<R>& original)
        : _bounds(original._bounds)
      { }
    
      /*! \brief Copy assignment operator. */
      Rectangle<R>& operator=(const Rectangle<R>& A) {
        if(this != &A) {
          this->_bounds = A._bounds;
        }
        return *this;
      }

      /*! \brief Assign from a rectangle expression. */
      template<class E>
      Rectangle<R>& operator=(const RectangleExpression<E>& original)
      {         
        const E& expression=original();
        this->_bounds.resize(2*expression.dimension());
        for (size_type i=0; i!=this->dimension(); ++i) {
          this->set_lower_bound(i,expression.lower_bound(i));
          this->set_upper_bound(i,expression.upper_bound(i));
        }
        return *this;
      }
    
      //@}
      
      
      //@{
      //! \name Conversion operators
      /*! \brief Convert to a Parma Polyhedra Library closed polyhedron. */
      operator Parma_Polyhedra_Library::C_Polyhedron () const;
      //@}
      
      
      //{@
      //! \name Comparison operators
      /*! \brief The equality operator */
      bool operator==(const Rectangle<R>& A) const
      {
        if (A.empty() && this->empty()) { return true; }
        if (A.empty() || this->empty()) { return false; }
        if(this->dimension()!=A.dimension()) { return false; }
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          if (this->lower_bound(i)!=A.lower_bound(i)) { return false; }
          if (this->upper_bound(i)!=A.upper_bound(i)) { return false; }
        }
        return true;
      }
      
      /*! \brief The inequality operator */
      bool operator!=(const Rectangle<R>& A) const {
        return !(*this == A);
      }
      //@}
      

      //@{
      //! \name Data access
      /*! \brief Returns the projection onto the \a i th coordinate. */
      Interval<R>& operator[] (dimension_type i) {
        return reinterpret_cast<Interval<R>&>(this->_bounds[2*i]);
      }
      //IntervalReference<R> operator[] (dimension_type i) {
      //  return IntervalReference<R>(this->_bounds[i]);
      //}
      
      /*! \brief Returns the projection onto the \a i th coordinate. */
      const Interval<R>& operator[] (dimension_type i) const {
        return reinterpret_cast<const Interval<R>&>(this->_bounds[2*i]);
      }
      
      /*! \brief Returns the lower bound of the \a i th coordinate */
      const Interval<R>& interval(dimension_type i) const {
        return reinterpret_cast<const Interval<R>&>(this->_bounds[2*i]);
      }
      
      /*! \brief Returns the lower bound of the \a i th coordinate */
      const R& lower_bound(dimension_type i) const {
        return this->_bounds[2*i];
      }
      
      /*! \brief Returns the upper bound of the \a i th coordinate */
      const R& upper_bound(dimension_type i) const {
        return this->_bounds[2*i+1];
      }
      
      /*! \brief The lower corner. */
      Point<R> lower_corner() const {
        Point<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->lower_bound(i);
        }
        return result;
      }
      
      /*! \brief The upper corner. */
      Point<R> upper_corner() const {
        Point<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->upper_bound(i);
        }
        return result;
      }
      
      /*! \brief The set of position vectors of the rectangle. */
      LinearAlgebra::Vector< Interval<R> > position_vectors() const {
        LinearAlgebra::Vector< Interval<R> > result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->interval(i);
        }
        return result;
      }
      //@}
      
      
      //@{ 
      //! \name Modifying operations
      /*! \brief Makes the rectangle empty. */
      void clear() {
        if(this->_bounds.size()!=0) {
          this->_bounds[0]=1;
          this->_bounds[1]=0;
        }
      }
      
      /*! \brief Sets the \a i th interval. */
      void set_interval(dimension_type i, Interval<R> x) {
        this->set_lower_bound(i,x.lower());
        this->set_upper_bound(i,x.upper());
      }
      
      /*! \brief Sets the lower bound of the \a i th coordinate to \a r. */
      void set_lower_bound(dimension_type i, const R& l) {
        this->_bounds[2*i]=l;
      }
      
      /*! \brief Sets the upper bound of the \a i th coordinate to \a u. */
      void set_upper_bound(dimension_type i, const R& u) {
        this->_bounds[2*i+1]=u;
      }

      /*! \brief Expand the Rectangle by \a delta in each direction. */
      Rectangle<R>& expand_by(const R& delta);
      //@}
      
      
      //@{
      //! \name Rectangle geometric operations
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const {
        return this->_bounds.size()/2;
      }
      
      /*! \brief True if the rectangle is empty. A zero-dimensional rectangle is considered empty. */
      bool empty() const {
        if(this->dimension()==0) {
          return true;
        }
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          if(this->lower_bound(i) > this->upper_bound(i)) {
            return true;
          }
        }
        return false;
      }
      
      /*! \brief True if the rectangle has empty interior. */
      bool empty_interior() const {
        if(this->dimension()==0) {
          return true;
        }
        for(size_type i=0; i!=this->dimension(); ++i) {
          if(this->lower_bound(i) >= this->upper_bound(i)) {
            return true;
          }
        }
        return false;
      }
      
      /*! \brief The centre. */
      Point<R> centre() const {
        Point<R> result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=Numeric::div_approx(Numeric::add_approx(this->lower_bound(i),this->upper_bound(i)),R(2));
        }
        return result;
      }
      
      /*! \brief The radius in the sup norm. */
      R radius() const {
        R diameter=0;
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          diameter=Numeric::max_up(diameter,Numeric::sub_up(this->upper_bound(i),this->lower_bound(i)));
        }
        return div_up(diameter,R(2));
      }
      
      /*! \brief Compute a quadrant of the Rectangle determined by \a q.
       *  \a q is a binary word such that the ith bit of q is 0 if the lower half
       *  of the rectangle in the ith coordinate is used, and 1 if the upper
       *  half is used.
       */
      Rectangle<R> quadrant(const Combinatoric::BinaryWord& q) const;
      /*! \brief Subdivide into smaller pieces. */
      ListSet<R,Geometry::Rectangle> subdivide() const;
      
      /*! \brief The vertices of the rectangle. */
      PointList<R> vertices() const;

      /*! \brief The \a i th vertex. */
      Point<R> vertex(size_type i) const;
        
      /*! \brief Tests if \a point is included into a rectangle. */
      bool contains(const Point<R>& pt) const;
      /*! \brief Tests if \a point is included into the interior a rectangle. */
      bool interior_contains(const Point<R>& pt) const;
      
      /*! \brief A rectangle containing the given rectangle; returns a copy. */
      Rectangle bounding_box() const {
        return *this;
      }
      //@}
      
#ifdef DOXYGEN
      //@{ 
      //! \name Binary geometric predicates
      /*! \brief Set equality operator. */
      friend bool equal(const Rectangle<R>& A, const Rectangle<R>& B) const;
     
      /*! \brief Tests disjointness with \a r. */
      friend bool disjoint(const Rectangle<R>& A, const Rectangle<R>& B) const;
      /*! \brief Tests intersection of interiors of \a A and \a B. */
      friend bool interiors_intersect(const Rectangle<R>& A, const Rectangle<R>& B) const;
      /*! \brief Tests inclusion in the interior of \a r. */
      friend bool inner_subset(const Rectangle<R>& A, const Rectangle<R>& B) const;
      /*! \brief Tests if the rectangle is a subset of another rectangle \a r. */
      friend bool subset(const Rectangle<R>& A, const Rectangle<R>& B) const;
      //@{ 
      //! \name Binary geometric operations
      /*! \brief The intersection of \a A and \a B. */
      friend Rectangle<R> intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
      friend Rectangle<R> regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The smallest rectangle containing \a A and \a B. */
      friend Rectangle<R> rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise sum of rectangles \a A and \a B. */
      friend Rectangle<R> minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise difference of rectangles \a A and \a B. */
      friend Rectangle<R> minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B); 
      
      /*! \brief The difference between two rectangles. */
      friend LinearAlgebra::Vector< Interval<R> > operator-(const Rectangle<R>& A, 
                                                        const Rectangle& B);
      /*! \brief Adds a vector to a rectangle. */
      friend Rectangle<R> operator+(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector<R>& v);
      /*! \brief Adds an interval vector to a rectangle. */
      friend Rectangle<R> operator+(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector< Interval<R> >& v);
      /*! \brief Subtracts a vector from a rectangle. */
      friend Rectangle<R> operator-(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector<R>& v);
      /*! \brief Subtracts an interval vector from a rectangle. */
      friend Rectangle<R> operator-(const Rectangle<R>& r, 
                                    const LinearAlgebra::Vector< Interval<R> >& v);
      //@}
#endif
      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     private:
      array<R> _bounds;
    };
  

    
    template <typename R>
    inline
    bool 
    Rectangle<R>::contains(const Point<R>& p) const 
    {
      const Rectangle<R>& self=*this;
      if (p.dimension()!=self.dimension()) {
        throw std::domain_error("The point has a different dimension to the rectangle.");
      }
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(!self[i].contains(p[i])) {
          return false;
        }
      }
      return true;
    }

    template <typename R>
    inline
    bool 
    Rectangle<R>::interior_contains(const Point<R>& p) const 
    {
      const Rectangle<R>& self=*this;
      if (p.dimension()!=self.dimension()) {
        throw std::domain_error("The point has a different dimension to the rectangle.");
      }
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(self.lower_bound(i) >= p[i] || self.upper_bound(i) <= p[i]) {
          return false;
        }
      }
      return true;
    }
      
    template <typename R>
    inline
    bool 
    equal(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if(A.dimension()!=B.dimension()) {
        return false;
      }
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(!Numeric::equal(A[i],B[i])) {
          return false;
        }
      }
      return true;
    }
      
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if(A.dimension()!=B.dimension()) {
        throw std::domain_error("disjoint(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(Numeric::disjoint(A[i],B[i])) {
          return true;
        }
      }
      return false;
    }
   
    template <typename R>
    inline
    bool 
    interiors_intersect(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("interiors_intersect(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(!(A[i].upper()>B[i].lower() && A[i].lower()<B[i].upper())) {
          return false;
        }
      }
      return true;
    }
    
    template <typename R>
    inline
    bool 
    inner_subset(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("inner_subset(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for (size_type i=0; i!=A.dimension(); ++i) {
        if(!(A[i].lower()>B[i].lower() && A[i].upper()<B[i].upper())) {
          return false;
        }
      }
      return true;
    }
    
    template <typename R>
    inline
    bool 
    subset(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("subset(Rectangle,Rectangle): The two rectangles have different dimensions");
      }
      for (size_type i=0; i< A.dimension(); ++i) {
        if(!Numeric::subset(A[i],B[i])) {
          return false;
        }
      }
      return true;
    }
    
    
    template <typename R>
    inline
    Rectangle<R> 
    intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("Cannot intersect with a rectangle of a different dimension.");
      }
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::intersection(A[i],B[i]);
      }
      return C;
    }

    template <typename R>
    inline
    Rectangle<R> 
    regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("Cannot intersect with a rectangle of a different dimension.");
      }
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::intersection(A[i],B[i]);
        if(C[i].lower()>=C[i].upper()) {
          C[i]=Interval<R>();
        }
      }
      return C;
    }

    template <typename R>
    inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error("The two parameters have different space dimensions");
      }      
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::hull(A[i],B[i]);
      }
      return C;
    }

    template<typename R> 
    Rectangle<R> 
    minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error(
              "minkowski_sum: the two rectangles have different dimension.");
      }
      for(dimension_type i=0; i!=C.dimension(); ++i) {
        C[i]=A[i]+B[i];
      }
      return C;
    }

    template<typename R> 
    Rectangle<R> 
    minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error(
              "minkowski_difference: the two rectangles have different dimension.");
      }
      for(dimension_type i=0; i!=C.dimension(); ++i) {
          C[i]=A[i]-B[i];
      }
      return C;
    }
    
    
    template <typename R> 
    inline
    bool 
    subset(const Rectangle<R>& A, ListSet<R,Geometry::Rectangle>& B);
    

    template <typename R> 
    inline
    bool 
    inner_subset(const Rectangle<R>& A, 
                 const ListSet<R,Geometry::Rectangle>& B);
    

    template <typename R> 
    inline
    bool 
    subset_of_open_cover(const Rectangle<R>& A, 
                         const ListSet<R,Geometry::Rectangle>& U);
    
    
    
    template<typename R>
    inline
    LinearAlgebra::Vector< Interval<R> > 
    operator-(const Geometry::Rectangle<R>& r1, 
              const Geometry::Rectangle<R>& r2)
    {
      assert(r1.dimension()==r2.dimension());
      
      return r1.position_vectors()-r2.position_vectors();
    }

    template<typename R, typename E>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const boost::numeric::ublas::vector_expression<E>& v)
    {
      const E& ev=v();
      LinearAlgebra::Vector< Interval<R> > iv=ev;
      return r+iv; 
    }
      
    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }

    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }

    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }

    template<typename R>
    inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==v.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }

    template<typename R>
    inline
    Geometry::Rectangle<R> 
    scale(const Geometry::Rectangle<R>& r, const R& scale_factor) {

      Geometry::Rectangle<R> result(r.dimension());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,scale_factor*r[i]);
      }

      return result;
    }
    
    template<typename R> inline 
    std::ostream& operator<<(std::ostream& os, const Rectangle<R>& r) {
      return r.write(os);
    }
    
    template<typename R> inline
    std::istream& operator>>(std::istream& is, Rectangle<R>& r) {
      return r.read(is);
    }

    
  }
}

#endif /* _ARIADNE_RECTANGLE_H */
