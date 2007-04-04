/***************************************************************************
 *            test_polyhedron.cc
 *
 *  Copyright  2005-6  Alberto Casagrande,  Pieter Collins
 *  Email casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>

#include <vector>

#include "test_float.h"

#include "ariadne.h"
#include "numeric/rational.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "output/epsfstream.h"

#include "geometry/ppl_polyhedron.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_polyhedron();
template<> int test_polyhedron<Rational>();

int main() {
  
  cout << boolalpha;
  
  test_polyhedron<Float>();
  test_polyhedron<Rational>();
   
  cerr << "INCOMPLETE ";
}



template<class R>
int 
test_polyhedron() 
{
  cout << "test_polyhedron<" << name<R>() << ">" << endl;
  Matrix<R> A("[1.0,0.875;-1,1.125;0.125,-2.25]");
  Vector<R> b("[1.375,0.5,0.25]");
  Polyhedron<R> phd1;
  
  phd1=Polyhedron<R>(A,b);
  cout << "phd1=" << phd1 << endl;
  cout << "phd1.dimension()=" << phd1.dimension() << ", " 
       << "phd1.number_of_constraints()=" << phd1.number_of_constraints() << endl;
  //assert(!phd1.empty());

  Point<R> pt1("(0.25,0.375)");
  cout << "pt1=" << pt1 << endl;

  typename Polyhedron<R>::constraints_const_iterator iter=phd1.constraints_begin();
  for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd1.constraints_begin();
      c_iter!=phd1.constraints_end(); ++c_iter)
  {
    const Constraint<R>& c=*c_iter;
    cout << c << flush;
    cout << "  " << c.satisfied_by(pt1) << endl;
  }
  
  cout << "phd1.contains(pt1)=" << flush; cout << phd1.contains(pt1) << endl;
  assert((bool)(phd1.contains(pt1)));
  
  Point<R> pt2("(2.25,-0.375)");
  cout << "pt2=" << pt2 << endl;
  for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd1.constraints_begin();
      c_iter!=phd1.constraints_end(); ++c_iter)
  {
    const Constraint<R>& c=*c_iter;
    cout << c << flush;
    cout << "  " << c.satisfied_by(pt2) << endl;
  }
  cout << "phd1.contains(pt2)=" << flush; cout << phd1.contains(pt2) << endl;
  assert((bool)(!phd1.contains(pt2)));
  

  Rectangle<R> r1("[-0.06125,0.25]x[0.125,0.375]");
  cout << "r1=" << r1 << endl;
  for(class Rectangle<R>::vertices_const_iterator v_iter=r1.vertices_begin();
      v_iter!=r1.vertices_end(); ++v_iter)
  {
    for(typename Polyhedron<R>::constraints_const_iterator c_iter=phd1.constraints_begin();
        c_iter!=phd1.constraints_end(); ++c_iter)
    {
      cout << *c_iter << ".satisfied_by" << *v_iter << "=" << c_iter->satisfied_by(*v_iter) << endl;
    }
  }
  cout << "subset(r1,phd1)=" << subset(r1,phd1) << endl;
  assert((bool)(subset(r1,phd1)));
  
  Rectangle<R> r2("[-0.125,0.25]x[0.125,0.75]");
  cout << "r2=" << r2 << endl;
  cout << "subset(r2,phd1)=" << subset(r2,phd1) << endl;
  assert((bool)(!subset(r2,phd1)));
  
  Zonotope<R> z1(r1);
  cout << "z1=" << z1 << endl;
  cout << "subset(z1,phd1)=" << subset(z1,phd1) << endl;
  //assert(subset(z1,phd1));
  
  Polytope<R> p2(r1);
  cout << "p2=" << flush; cout << p2 << endl;
  cout << "subset(p2,phd1)=" << subset(p2,phd1) << endl;
  assert((bool)(subset(p2,phd1)));
  
  cout << endl;


  Polyhedron<R> phd2(Matrix<R>("[1,-1;1,-2;-1,1;-1,2]"),Vector<R>("[1.375,1.625,0.625,0.375]"));
  cout << "phd2=" << phd2 << endl;
  cout << "phd2.bounding_box()=" << phd2.bounding_box() << endl;

  Polytope<Rational> qpltp2=Polytope<Rational>(Polyhedron<Rational>(phd2));
  cout << "qpltp2=" << qpltp2 << endl << "qpltp2.bounding_box()=" << qpltp2.bounding_box() << endl;
  Matrix<R>  G=Matrix<R>("[0.889,0.615,-1.600,-2.667;-0.889,-1.231,1.600,5.333]");
  G=Matrix<R>("[3.125,1.125,-2.875,-0.875;1.75,-0.25,-2.25,-0.25]");
  cout << "G=" << G << endl;
  Polytope<R> pltp2=Polytope<R>(G);
  cout << "pltp2=" << pltp2 << endl << "pltp2.bounding_box()=" << pltp2.bounding_box() << endl;

  Rectangle<R> bbox2=phd2.bounding_box().expand_by(0.25);
  bbox2=Rectangle<R>("[-4,4]x[-4,4]");
  Grid<R> gr2(2,0.125);
  GridCellListSet<R> uap2=under_approximation(phd2,gr2);
  GridCellListSet<R> oap2=over_approximation(phd2,gr2);
  cout << "uap2.size()=" << uap2.size() << endl;
  cout << "oap2.size()=" << oap2.size() << endl;
  epsfstream eps;
  eps.open("test_polyhedron-1.eps",bbox2);
  eps.set_fill_colour("white");
  eps << bbox2;
  eps.set_fill_colour("red");
  eps << oap2;
  eps.set_fill_colour("green");
  eps << phd2;
  eps.set_fill_colour("blue");
  eps << uap2;
  eps.close();
  
  cout << endl;
  
  return 0;
}


template<>
int 
test_polyhedron<Rational>() 
{
  typedef Rational R;
  cout << "test_polyhedron<" << name<R>() << ">" << endl;
  LinearAlgebra::Matrix<R> A("[1,7/8;-1,9/8;1/8,-9/4]");
  LinearAlgebra::Vector<R> b("[11/8,1/2,1/4]");
  Polyhedron<R> phd(A,b);
  cout << "phd=" << phd << endl;
  cout << "phd.constraints()=" << phd.constraints() << endl;
    
  Polytope<R> pltp(phd);
  cout << "Polytope(phd)=" << pltp << endl;
  cout << "  Vertices should be (71/128,15/16), (-2/3,-4/7),  (212/151,-5/151)" << endl;
  Polyhedron<R> phd2(pltp);
  cout << "Polyhedron(Polytope(phd))=" << phd2 << endl;
   cout << endl;
 
  Rectangle<R> r1("[7/8,8/8]x[6/8,7/8]");
  Rectangle<R> r2("[2/8,3/8]x[5/8,6/8]");
  Rectangle<R> r3("[3/8,4/8]x[4/8,5/8]");

  Rectangle<R> bbox=phd.bounding_box().expand_by(R(0.5));

  epsfstream eps;
  eps.open("test_polyhedron-2.eps",bbox);
  eps.set_fill_colour("white");
  eps << bbox;
  eps.set_fill_colour("green");
  eps << phd;
  eps.set_fill_colour("red");
  eps << r1;
  eps.set_fill_colour("yellow");
  eps << r2;
  eps.set_fill_colour("blue");
  eps << r3;
  eps.close();

  cout << "phd=" << phd << endl;
  cout << "r1=" << r1 << endl;
  cout << "r2=" << r2 << endl;
  cout << "r3=" << r3 << endl;
  cout << "disjoint(r1,phd)=" << disjoint(r1,phd) << endl;
  cout << "disjoint(r2,phd)=" << disjoint(r2,phd) << endl;
  cout << "subset(r3,phd)=" << subset(r3,phd) << endl;
  assert((bool)(disjoint(r1,phd)==true));
  assert((bool)(disjoint(r2,phd)==false));
  assert((bool)(subset(r3,phd)==true));

  return 0;
}
