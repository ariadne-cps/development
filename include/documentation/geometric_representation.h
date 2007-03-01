/***************************************************************************
 *            geometric_representation.h
 *
 *  Copyright  2004-7  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

/*! 

\file geometric_representation.h
\brief Documentation for geometric sets and operations



\page geometric Geometric Representation

\deprecated The material in this section is out-of-date

\section Introduction

One of the main concerns of Ariadne is the computation of points and sets in Euclidean space.
Unfortunately, the set of points in Euclidean space is uncountable (it has continuum cardinality) and so there
is no way of modelling all points in Euclidean space using binary words. The set of subsets of Euclidean space
is also uncountable, and has an even greater cardinality!

To represent an uncountable set, we can use approximations.
We take a set of elements which we model exactly, such as \a double reals or \a rational reals, which we call <em>denotable</em> elements.
This denotable set may be finite or countably infinite, the former being denoted by words of a fixed length, and the latter by arbitrary length words.
An <em>approximation</em> to an element is then a pair consisting of a denotable element and an <em>error</em>.
Ideally, we wish to be able to approximate to arbitrary precision;
this means that the set of denotable elements must be <em>dense</em> in the set of all elements.
We can then represent arbitrary elements be a <em>convergent sequence</em> of approximations.
In terms of a concrete representation on a computer, we can think of a such a sequence as a neverending <em>stream</em> of data.

The material in this section is heavily influenced by the book "Computational Analysis" by Klaus Weihrauch.

\section state Representation of points in Euclidean space

Points in Euclidean space can be represented by the templated class \a DenotablePoint<T>, where \a T is a numeric type such as double or rational. From this, we can define a PointApproximation class, which gives an approximation of a state, plus an error bound in terms of the sup norm or Euclidean norm.

\code
concept State
{
  typename real_type; 

  State(const State &);
  State& operator=(const State &);
 
  real_type operator[] (size_type) const;
};
\endcode 

\section basicset Basic sets.

A \c BasicSet provides a building block for representing more complicated sets.
Mathematically, a \c BasicSet type represents elements of a countable base of a topological space.
More precisely, a \c BasicSet represents the \em closure of a basic set for the topology.

In many cases, it is not possible to evaluate geometric predicates concerning two sets using a given real number type. For this reason, the results of a geometric predicate returns an object of type \a tribool, which may be \a true, \a false or \a indeterminate (unknown).
 


\code
// The basic set concept.
concept BasicSet
{
  type real_type; // The type of denotable real number used for the representation.
  type state_type; // The type of denotable point the set contains.

  BasicSet(const std::string &); // Construct from a string literal (optional).
  BasicSet(const BasicSet &); // Copy constructor.
  BasicSet & operator=(const BasicSet &); // Assignment operator.

  dimension_type dimension() const; // The dimension of the set.
  state_type centre() const; // A point in the set (typically, the "centre" point, if this makes sense).
  real_type radius() const; // The maximum distance from the centre to another point in the set in an appropriate metric. 
  real_type volume() const; // An approximation to the volume of the set. (Optional).

  tribool empty() const; // Tests if the set is empty.
  tribool bounded() const; // Tests if the set is bounded.

  tribool contains(const State &) const; // Tests if the set contains a point.

  Rectangle<real_type> bounding_box() const; // A rectangle containing the set.
  tribool disjoint(const Rectangle &); // Tests if the set is disjoint from a Rectangle
  tribool superset(const Rectangle &); // Tests if the set contains a Rectangle.

};

 tribool equal(const BasicSet1 &, const BasicSet2 &); // Returns indeterminate if the two sets are equal; if false, then they are not equal.
 tribool disjoint(const BasicSet1 &, const BasicSet2 &); // If true, then the sets are disjoint; if false, then they robustly intersect.
 tribool subset(const BasicSet1 &, const BasicSet2 &); // If true, then the first set is a subset of the second interior of the second; 
                                                     // if false, then the first set is not a subset of the second.

 // Optional, depending on whether the operation yields a basic set of the same type.
 BasicSet open_intersection(const BasicSet &, const BasicSet &); // The closure of the intersection of the interiors of the two sets.
 BasicSet closed_intersection(const BasicSet &, const BasicSet &); // The intersection of the two (closed) sets.
 BasicSet convex_hull(const BasicSet &, const BasicSet &); // The convex hull of the two sets.
 BasicSet minkowski_sum(const BasicSet &, const BasicSet &); // The Minkowski (pointwise) sum of the two sets.
 BasicSet minkowski_difference(const BasicSet &, const BasicSet &); // The Minkowski (pointwise) difference of the two sets.

\endcode

Classes fulfilling the BasicSet concept are \ref Rectangle (or Cuboid), Simplex, Parallelotope, Zonotope, Polytope, Polyhedron, Sphere and Ellipsoid.
Actually, these are templates, parameterised by the real number type real_type.


\section denotable_set Denotable Sets

A DenotableSet implements a set as a union of basic sets type \c DenotableSet::basic_set_type.
\code
concept DenotableSet
{
  type real_type;
  type state_type;
  type basic_set_type;

  type const_iterator; // Must satisfy the requirements of a ForwardIterator.

  // No default constructor required.

  DenotableSet(const DenotableSet &);
  DenotableSet & operator=(const DenotableSet &);

  // No equality operator required.

  // SetInterface-theoretic operations
  dimension_type dimension() const;
  tribool empty() const;
  tribool contains(const state_type &) const;

  Rectangle<real_type> bounding_box() const; // Optional.

  void adjoin(const basic_set_type &);
  void adjoin(const DenotableSet &);

  // List operations
  const_iterator begin() const;
  const_iterator end() const;

  size_type size() const; // Only required if the iterator is a RandomAccessIterator.
  basic_set_type operator[] (size_type) const; // Only required if the iterator is a RandomAccessIterator.

  void push_back(const basic_set_type &); // Only used if the DenotableSet is an ordered list. (Optional)
  basic_set_type pop_back(); // Only used if the DenotableSet is an ordered list. (Optional)

  void insert(const basic_set_type &); // Only used if the DenotableSet is an unordered or sorted list. (Optional)
  void remove(const basic_set_type &); // Only used if the DenotableSet is an unordered or sorted list. (Optional)
};

tribool subset(const BasicSet &, const DenotableSet &); // Optional, but highly recommended.

tribool disjoint(const DenotableSet &, const DenotableSet &);
tribool subset(const DenotableSet &, const DenotableSet &);

DenotableSet join(const DenotableSet &, const DenotableSet &);
DenotableSet open_intersection(const BasicSet &, const DenotableSet &); // Optional.
DenotableSet difference(const DenotableSet&, const DenotableSet&); // Optional

\endcode

\section set_approximation Approximating Sets

%Ariadne provides operators for approximating sets. All the operators have one
of the following forms.

\code
Result outer_approximation(Argument,Error); // postcondition: inner_subset(Argument,Result)
Result over_approximation(Argument,Error);  // postcondition: subset(Argument,Result)
Result lower_approximation(Argument,Error); // postcondition: 
Result under_approximation(Argument,Error); // postcondition: subset(Result,Argument)
Result inner_approximation(Argument,Error); // postcondition: inner_subset(Result,Argument)
\endcode

The error specification depends on the type of approximation used. 
 - When approximating by a GridSet, the error is determined by the grid.
 - When approximating by a PartitionTreeSet, the error is determined by the 
      partition scheme and the depth.
 - When approximating by a ListSet, the error is a metric error bound.








\section GeometricOps Geometric Operations

The core geometric types used by %Ariadne to represent are Rectangle, Zonotope, Polytope (described by generators)
and Polyhedron (described by constraints). 
The core geometric operations are contains(A,p), subset(A,B) and disjoint(A,B) (equivalent to intersects(A,B) ) .

  - Rectangle representation: \f$l\leq x\leq u\f$.
  - Zonotope representation: \f$x=c+Ge,\ -1\leq e\leq1\f$.
  - Polytope representation: \f$x=Vs,\ \sum s_i=1,\ s\geq0\f$.
  - Polyhedron representation: \f$Ax\leq b\f$.


Alternative representations for polytopes and polyhedron are given by taking \f$\hat{x}=(x,1)\f$.
  - Rectangle representation: \f$x_i\in[l_i,u_i]\f$.
  - Zonotope representation: \f$x = \hat{G}\hat{e},\ e_0=1,\ -1\leq e\leq1\f$.
  - Polytope representation: \f$\hat{x} = Rs\f$, \f$s\geq0\f$.
  - Polyhedron representation: \f$C\hat{x}\geq0\f$.

We sometimes use the based polytopic representation
  - Polytope representation: \f$x=v+Ws,\ \sum s_i\leq1,\ s\geq0\f$.

Henceforth, we shall always assume \f$-1\leq e\leq1\f$ and \f$1\!\cdot\!s=1;\ s\geq0\f$.


\section geometricpreprocessing Preprocessing sets

To simplify certain computations, it may be useful to pre-process the matrices describing a Zonotope, Polytope or Polyhedron.

Given an \f$m\times n\f$ matrix \f$A\f$ of rank \f$k\f$, find an index set \f$I\subset\{0,1,\ldots,n\!-\!1\}\f$ of cardinality \f$k\f$, and an \f$k\times m\f$ matrix \f$B^{-1}\f$ such that
\f[ B^{-1} A_I = \mathbf{I} \f]
where \f$ A_I \f$ is the \f$m\times k\f$ matrix consisting of the columns of \f$A\f$ with index \f$j\in I\f$.


\section conversion Converting between basic sets

A %Rectangle can be directly converted to a %Polytope or %Polyhedron without using arithmetic. 
A %Rectangle can be easily converted to a %Zonotope, and a %Zonotope to a %Polytope, but these conversions require arithmetic.
Conversion between a %Polytope to a %Polyhedron can be performed using the double description algorithm.

The conversion from a %Rectangle or %Zonotope to a %Polytope or conversion between %Polytope and %Polyhedron may be of exponential complexity.

 \section contains Testing inclusion

  - contains(Rectangle,Point) : Check \f$l\leq x\leq u\f$.
  - contains(Zonotope,Point) : Solve \f$x=c+Ge;\ -1\leq e\leq1\f$.
  - contains(Polytope,Point) : Solve \f$x=Vs;\ 1\!\cdot\!s=1;\ s\geq 0\f$.
  - contains(Polyhedron,Point) : Check \f$Ax\leq b\f$.

 \section Intersection Testing intersection/disjointness

  - intersects(Rectangle,Rectangle) : Check \f$l_1\leq u_2;\ l_2\leq u_1\f$.
  - intersects(Rectangle,Zonotope) : Solve \f$x=c+Ge;\ l\leq x\leq u;\ -1\leq e\leq1\f$ or \f$l\leq c+Ge\leq u;\ -1\leq e\leq 1\f$.
  - intersects(Zonotope,Zonotope) : Solve \f$c_1+G_1e_1 = c_2+G_2e_2;\ -1\leq e_1\leq1;\ -1\leq e_2\leq1\f$.
  - intersects(Rectangle,Polytope) : Solve \f$x=Vs;\ 1\!\cdot\!s=1;\ l\leq x\leq u;\ s\geq 0\f$ or \f$l\leq Vs\leq u;\ 1\!\cdot\!s=1;\ s\leq 1\f$.
  - intersects(Zonotope,Polytope) : Solve \f$c+Ge = Vx;\ \ 1\!\cdot\!s=1;\ -1\leq e\leq 1;\ s\geq0\f$.
  - intersects(Polytope,Polytope) : Solve \f$V_1s_1=V_2s_2;\ 1\!\cdot\!s_1=1;\ 1\!\cdot\!s_2=1;\ \cdot s_1\geq0;\ s_2\geq0\f$.

  - intersects(Rectangle,Polyhedron) : Solve \f$Ax\leq b;\ l\leq x\leq u\f$.
  - intersects(Zonotope,Polyhedron) : Solve \f$Ac+AGe \leq b;\ -1\leq e\leq 1\f$.
  - intersects(Polytope,Polyhedron) : Solve \f$AVs\leq b;\ 1\!\cdot\!s=1;\ s\geq0\f$.
  - intersects(Polyhedron,Polyhedron) : Solve \f$A_1x\leq b_1;\ A_2x\leq b_2\f$.

All these problems can be solved using (a variant of) the simplex algorithm. 
When computing over-approximations of a set on a grid using disjoint(Rectangle R, ConvexSet S), it may be useful to pre-process S to simplify computation of an initial vertex for the feasibility problem.


\section Subset Testing subset

Testing whether a convex set is a subset of a polyhedron, or whether a polytope is a subset of a convex set are easy:
  - subset(ConvexSet S, Polyhedron P) : Check disjoint(S,H) for all complementary halfspaces H of P.
  - subset(Polytope P, ConvexSet S) : Check contains(S,v) for all vertices v of P.

Checking whether a Rectangle, Zonotope or Polytope is a subset of a Rectangle can be performed directly:
  - subset(Rectangle, Rectangle) : Check \f$l_1\geq l_2\ \textrm{and}\ u_1\leq u_2\f$.
  - subset(Zonotope, Rectangle) : Check \f$\sum_{j} |G_{ij}|\leq \min\{c_i-l_i,u_i-c_i\} \forall\,i\f$.
  - subset(Polytope, Rectangle) : Check \f$l\leq v\leq u\f$ for all vertices \f$v\f$ of \f$P\f$.

The following operations cannot be performed efficiently, but can be easily implemented by testing the vertices of the Rectangle, and are provided for conformance to the SetInterface interface.
  - subset(Rectangle, Zonotope)
  - subset(Rectangle, Polytope)

The following operation is provided to help with checking zonotopic over-approximations.
  - subset(Zonotope, Zonotope)

While a Rectangle can be directly converted to both a Polytope and a Polyhedron, the conversion to a Polytope yields \f$2^d\f$ vertices.
Hence testing subset(Rectangle, ConvexSet) is in general of exponential complexity, as is subset(Zonotope, ConvexSet).
When testing subset(Rectangle R, ConvexSet S) for under-approximationg S on a grid, we need to check contains(S,v) for each grid vertex v once only, considerably reducing the amount of work.

Note that for a single test of subset(Rectangle R, Polyhedron P), we prefer to use disjoint(R,H), but for under-approximating a polyhedron P on a grid, we prefer to check contains(P,v) for the grid vertices.


\section PolyhedralConversion Converting between a polyhedron and a polytope.
 
Augment the state by \f$\hat{x}=(x,1)\f$.

  - Define a Polytope by generators \f$x=\lambda g\f$, given as columns of the augmented generator matrix \f$G'\f$.
  - Define a Polyhedron by constraints \f$a\cdot x\geq 0\f$, given as rows of the augmented constraints matric \f$A'\f$.

Construct the <em>saturation matrix</em> by 
  - Generator \f$g\f$ \e violates constraint \f$a\f$ if \f$a\cdot g<0\f$, 
  - Generator \f$g\f$ \e saturates constraint \f$a\f$ if \f$a\cdot g=0\f$.
  - Generator \f$g\f$ \e satisfies constraint \f$a\f$ if \f$a\cdot g=0\f$.
Saturation matrix \f$S=\mathrm{sgn}(AG)\f$.

Generators are \em adjacent if the corresponding columns of the saturation matrix
differ only in one row.

\section zonotope Zonotopic reduction methods

Throughout this sections, we use the supremum norm on \f$R^n\f$, and the corresponding operator norm on \f$\mathbb{R}^{m\times n}\f$.

Given a zonotope \f$ Z=\{ c+Ae \mid ||e||\leq 1 \}\subset \mathbb{R}^n\f$, where \f$A\in \mathbb{R}^{n\times p}\f$, 
we wish to compute a zonotope \f$Z' = \{ c + A' e' \mid ||e'||\leq 1\}\f$ with fewer generators 
i.e. \f$A'\in \mathbb{R}^{n\times p'}\f$ with \f$p'<p\f$.
The general reduction method is to choose \f$ A'\f$ such that \f$ A = A' B\f$ with \f$||B||\leq 1\f$.
The key to zonotopic reduction is to choose a method with good properties.

A simple criterion to note is that if \f$\sum_{j=1}^{p} |b_{ij}|<1\f$ for some \f$i\f$, then we can improve the approximation by taking
\f$D=\mathrm{diag}(d_{i})\f$ with \f$d_{i}=\sum_{j=1}^{p} |b_{ij}|\f$, 
and \f$B'= D^{-1}B\f$ which has \f$b'_{ij}=b_{ij}/\sum_{k=1}^{p} |b_{ik}|\f$.

It is clear that if the rows of \f$B\f$ are close to a set of mutually orthogonal coordinate vectors, then the approximation is good, 
since the image of \f$B\f$ is close to the unit ball. 


\subsection interval_zonotope Interval zonotopic reduction

An <em>interval zonotope</em> is a set of the form \f$ \{ y = c + A e \mid c\in R,\ A\in\mathcal{A} \text{ and } ||e||\leq1 \} 
= R + \mathcal{A} B\f$.
To reduce an interval zonotope, we first write \f$ R = \{ c + Be\mid ||e||\leq 1 \}\f$ and combine this in \f$\mathcal{A}\f$.
To reduce \f$ \mathcal{A} \f$, write \f$ \mathcal{A} = A \mathcal{B} \mathcal{A} \f$ where \f$ A\in\mathcal{A}\f$ and \f$A\mathcal{B}\ni I\f$.
Then take \f$ \mathcal{C} = \mathcal{B} \mathcal{A} \f$ and \f$ || \mathcal{C} || 
  = \sup_{i]1}^{n} \sum_{j=1}^{p} \max |\mathcal{C}_{ij}| \f$, 
where \f$ \max |\mathcal{C}_{ij}| = \max\{ |x| \mid x\in \mathcal{C}_{ij}\f$.

*/
