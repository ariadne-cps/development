/***************************************************************************
 *            geometry_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
 *
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

#include "config.h"

#include "geometry.h"
#include "point.h"
#include "curve.h"
#include "box.h"
#include "zonotope.h"
#include "polytope.h"
#include "polyhedron.h"
#include "grid_set.h"
#include "function_set.h"
#include "taylor_set.h"
#include "affine_set.h"

#include "discrete_location.h"
#include "hybrid_set.h"

#include "utilities.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {


template<>
struct from_python<Point> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Point>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        extract<boost::python::tuple> xtup(obj_ptr);
        extract<boost::python::list> xlst(obj_ptr);
        Point pt;
        if(xtup.check()) {
            boost::python::tuple tup=xtup(); pt=Point(len(tup));
            for(int i=0; i!=len(tup); ++i) { pt[i]=extract<double>(tup[i]); }
        } else if(xlst.check()) {
            boost::python::list lst=xlst(); pt=Point(len(lst));
            for(int i=0; i!=len(lst); ++i) { pt[i]=extract<double>(lst[i]); }
        }
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Point(pt);
        data->convertible = storage;
    }
};

template<>
struct from_python<Box> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Box>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        boost::python::list lst=extract<boost::python::list>(obj_ptr);
        Box* bx_ptr = new (storage) Box(len(lst));
        for(int i=0; i!=len(lst); ++i) { (*bx_ptr)[i]=extract<Interval>(lst[i]); }
        data->convertible = storage;
    }
};


template<class ES>
struct to_python< ListSet<ES> > {
    to_python() { boost::python::to_python_converter< ListSet<ES>, to_python< ListSet<ES> > >(); }

    static PyObject* convert(const ListSet<ES>& ls) {
        boost::python::list result;
        for(typename ListSet<ES>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class ES>
struct to_python< ListSet< HybridBasicSet<ES> > > {
    typedef ListSet< HybridBasicSet<ES> > SetType;
    to_python() { boost::python::to_python_converter< SetType, to_python<SetType> >(); }

    static PyObject* convert(const SetType& hls) {
        boost::python::dict result;
        for(typename SetType::locations_const_iterator iter=hls.locations_begin(); iter!=hls.locations_end(); ++iter) {
            result[iter->first]=iter->second;
        }
        return boost::python::incref(boost::python::dict(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};



class OpenSetWrapper
  : public OpenSetInterface, public wrapper< OpenSetInterface >
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool covers(const Box& r) const { return this->get_override("covers")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class ClosedSetWrapper
  : public ClosedSetInterface, public wrapper< ClosedSetInterface >
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool disjoint(const Box& r) const { return this->get_override("disjoint")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class OvertSetWrapper
  : public OvertSetInterface, public wrapper< OvertSetInterface >
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class CompactSetWrapper
  : public CompactSetInterface, public wrapper< CompactSetInterface >
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool disjoint(const Box& r) const { return this->get_override("disjoint")(); }
    tribool inside(const Box& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    Box bounding_box() const { return this->get_override("bounding_box")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class LocatedSetWrapper
  : public LocatedSetInterface, public wrapper< LocatedSetInterface >
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool covers(const Box& r) const { return this->get_override("covers")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    tribool disjoint(const Box& r) const { return this->get_override("disjoint")(); }
    tribool inside(const Box& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    Box bounding_box() const { return this->get_override("bounding_box")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

}


void export_set_interface() {
    class_<OpenSetInterface, boost::noncopyable> open_set_wrapper_class("OpenSetInterface", no_init);
    open_set_wrapper_class.def("covers",&OpenSetInterface::covers);
    open_set_wrapper_class.def("overlaps",&OpenSetInterface::overlaps);

    class_<CompactSetInterface, boost::noncopyable> compact_set_wrapper_class("CompactSetInterface", no_init);
    compact_set_wrapper_class.def("disjoint",&CompactSetInterface::disjoint);
    compact_set_wrapper_class.def("inside",&CompactSetInterface::inside);
    compact_set_wrapper_class.def("bounding_box",&CompactSetInterface::bounding_box);

    class_<LocatedSetInterface, boost::noncopyable> located_set_wrapper_class("LocatedSetInterface", no_init);

    class_<DrawableInterface,boost::noncopyable>("DrawableInterface",no_init);

}


void export_point()
{
    class_<Point> point_class("Point",init<Point>());
    point_class.def(init<uint>());
    point_class.def("__getitem__", &__getitem__<Point,int,Float>);
    point_class.def(self_ns::str(self));

    from_python<Point>();
    implicitly_convertible<Vector<Float>,Point>();

}

void export_box()
{
    typedef Vector<Interval> IVector;

    class_<Box,bases<CompactSetInterface,OpenSetInterface,Vector<Interval>,DrawableInterface > > box_class("Box",init<Box>());
    box_class.def(init<uint>());
    box_class.def(init< Vector<Interval> >());
    box_class.def("__eq__", (bool(*)(const Vector<Interval>&,const Vector<Interval>&)) &operator==);
    box_class.def("dimension", (uint(Box::*)()const) &Box::dimension);
    box_class.def("centre", (Point(Box::*)()const) &Box::centre);
    box_class.def("radius", (Float(Box::*)()const) &Box::radius);
    box_class.def("separated", (tribool(Box::*)(const Box&)const) &Box::disjoint);
    box_class.def("overlaps", (tribool(Box::*)(const Box&)const) &Box::overlaps);
    box_class.def("covers", (tribool(Box::*)(const Box&)const) &Box::covers);
    box_class.def("inside", (tribool(Box::*)(const Box&)const) &Box::inside);
    box_class.def("widen", (Box(Box::*)()const) &Box::widen);
    box_class.def("split", (std::pair<Box,Box>(Box::*)()const) &Box::split);
    box_class.def("split", (std::pair<Box,Box>(Box::*)(uint)const) &Box::split);
    box_class.def("split", (std::pair<Box,Box>(Box::*)()const) &Box::split);
    box_class.def(self_ns::str(self));

    def("split", (std::pair<Box,Box>(*)(const Box&)) &split);
    def("disjoint", (bool(*)(const IVector&,const IVector&)) &disjoint);
    def("subset", (bool(*)(const IVector&,const IVector&)) &subset);

    def("hull", (Box(*)(const Box&,const Box&)) &hull);
    def("intersection", (Box(*)(const Box&,const Box&)) &intersection);

    from_python<Box>();
    to_python< std::pair<Box,Box> >();
    implicitly_convertible<Vector<Interval>,Box>();

}

void export_zonotope()
{
    class_<Zonotope,bases<CompactSetInterface,OpenSetInterface,DrawableInterface> > zonotope_class("Zonotope",init<Zonotope>());
    zonotope_class.def(init< Point, Matrix<Float>, Vector<Float> >());
    zonotope_class.def(init< Point, Matrix<Float> >());
    zonotope_class.def(init< Box >());
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    def("contains", (tribool(*)(const Zonotope&,const Point&)) &contains);
    def("disjoint", (tribool(*)(const Zonotope&,const Box&)) &disjoint);
    def("overlaps", (tribool(*)(const Zonotope&,const Box&)) &overlaps);
    def("disjoint", (tribool(*)(const Zonotope&,const Zonotope&)) &disjoint);
}

void export_polytope()
{
    class_<Polytope,bases<LocatedSetInterface,DrawableInterface> > polytope_class("Polytope",init<Polytope>());
    polytope_class.def(init<int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));

}

void export_curve()
{
    to_python< std::pair<const Float,Point> >();

    class_<InterpolatedCurve,bases<DrawableInterface> > interpolated_curve_class("InterpolatedCurve",init<InterpolatedCurve>());
    interpolated_curve_class.def(init<Float,Point>());
    interpolated_curve_class.def("insert", &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__",boost::python::range(&InterpolatedCurve::begin,&InterpolatedCurve::end));
    interpolated_curve_class.def(self_ns::str(self));


}

void export_taylor_set()
{
    class_<TaylorImageSet,bases<CompactSetInterface,DrawableInterface> > taylor_set_class("TaylorImageSet",init<TaylorImageSet>());
    taylor_set_class.def(init<uint>());
    taylor_set_class.def(init<Box>());
    taylor_set_class.def(init< Vector<TaylorModel> >());
    //taylor_set_class.def(init<Zonotope>());
    taylor_set_class.def("domain", &TaylorImageSet::domain);
    taylor_set_class.def("models", &TaylorImageSet::models,return_value_policy<copy_const_reference>());
    taylor_set_class.def("bounding_box", &TaylorImageSet::bounding_box);
    taylor_set_class.def("range", &TaylorImageSet::bounding_box);
    taylor_set_class.def("split", (std::pair<TaylorImageSet,TaylorImageSet>(TaylorImageSet::*)()const) &TaylorImageSet::split);
    taylor_set_class.def("split", (std::pair<TaylorImageSet,TaylorImageSet>(TaylorImageSet::*)(uint)const) &TaylorImageSet::split);
    taylor_set_class.def(self_ns::str(self));

    def("split", (std::pair<TaylorImageSet,TaylorImageSet>(TaylorImageSet::*)()const) &TaylorImageSet::split);
    def("outer_approximation", (GridTreeSet(*)(const TaylorImageSet&,const Grid&,uint)) &outer_approximation);
    def("adjoin_outer_approximation", (void(*)(GridTreeSet&,const TaylorImageSet&,uint)) &adjoin_outer_approximation);
    def("zonotope", (Zonotope(*)(const TaylorImageSet&)) &zonotope);

    def("product", (TaylorImageSet(*)(const TaylorImageSet&,const Interval&)) &product);
    def("product", (TaylorImageSet(*)(const TaylorImageSet&,const Box&)) &product);
    def("product", (TaylorImageSet(*)(const TaylorImageSet&,const TaylorImageSet&)) &product);

    def("apply",(TaylorModel(*)(const ScalarFunction&,const TaylorImageSet&)) &apply);
    def("apply",(TaylorModel(*)(const ScalarTaylorFunction&,const TaylorImageSet&)) &apply);
    def("apply",(TaylorImageSet(*)(const VectorFunction&,const TaylorImageSet&)) &apply);
    def("apply",(TaylorImageSet(*)(const VectorTaylorFunction&,const TaylorImageSet&)) &apply);

    def("unchecked_apply",(TaylorImageSet(*)(const VectorTaylorFunction&,const TaylorImageSet&)) &unchecked_apply);

    implicitly_convertible<Box,TaylorImageSet>();

    to_python< std::pair<TaylorImageSet,TaylorImageSet> >();
}



void export_affine_set()
{

    class_<AffineSet,bases<DrawableInterface> >
        affine_set_class("AffineSet",init<AffineSet>());
    affine_set_class.def(init<Vector<Interval>, Matrix<Float>, Vector<Float> >());
    affine_set_class.def(init<Matrix<Float>, Vector<Float> >());
    affine_set_class.def("new_inequality_constraint", &AffineSet::new_inequality_constraint);
    affine_set_class.def("new_equality_constraint", &AffineSet::new_equality_constraint);
    affine_set_class.def("dimension", &AffineSet::dimension);
    affine_set_class.def("bounded", &AffineSet::bounded);
    affine_set_class.def("empty", &AffineSet::empty);
    affine_set_class.def("bounding_box", &AffineSet::bounding_box);
    affine_set_class.def("disjoint", &AffineSet::disjoint);
    affine_set_class.def("adjoin_outer_approximation_to", &AffineSet::adjoin_outer_approximation_to);
    affine_set_class.def("outer_approximation", &AffineSet::outer_approximation);
    affine_set_class.def("boundary", &AffineSet::boundary);
    affine_set_class.def(self_ns::str(self));
}


void export_constrained_image_set()
{
    class_<ConstrainedImageSet,bases<DrawableInterface> >
    constrained_image_set_class("ConstrainedImageSet",init<ConstrainedImageSet>());
    constrained_image_set_class.def(init<Box>());
    constrained_image_set_class.def(init<Box,VectorFunction>());
    constrained_image_set_class.def("domain", &ConstrainedImageSet::domain,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("function", &ConstrainedImageSet::function,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("constraint", &ConstrainedImageSet::constraint,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("number_of_parameters", &ConstrainedImageSet::number_of_parameters);
    constrained_image_set_class.def("number_of_constraints", &ConstrainedImageSet::number_of_constraints);
    constrained_image_set_class.def("apply", &ConstrainedImageSet::apply);
    constrained_image_set_class.def("new_space_constraint", (void(ConstrainedImageSet::*)(const NonlinearConstraint&))&ConstrainedImageSet::new_space_constraint);
    constrained_image_set_class.def("new_parameter_constraint", (void(ConstrainedImageSet::*)(const NonlinearConstraint&))&ConstrainedImageSet::new_parameter_constraint);
    //constrained_image_set_class.def("outer_approximation", &ConstrainedImageSet::outer_approximation);
    constrained_image_set_class.def("affine_approximation", &ConstrainedImageSet::affine_approximation);
    //constrained_image_set_class.def("affine_over_approximation", &ConstrainedImageSet::affine_over_approximation);
    constrained_image_set_class.def("adjoin_outer_approximation_to", &ConstrainedImageSet::adjoin_outer_approximation_to);
    constrained_image_set_class.def("inside", &ConstrainedImageSet::inside);
    constrained_image_set_class.def("disjoint", &ConstrainedImageSet::disjoint);
    constrained_image_set_class.def("overlaps", &ConstrainedImageSet::overlaps);
    constrained_image_set_class.def("split", (Pair<ConstrainedImageSet,ConstrainedImageSet>(ConstrainedImageSet::*)()const) &ConstrainedImageSet::split);
    constrained_image_set_class.def("split", (Pair<ConstrainedImageSet,ConstrainedImageSet>(ConstrainedImageSet::*)(uint)const) &ConstrainedImageSet::split);
    constrained_image_set_class.def(self_ns::str(self));

}


void export_hybrid_box()
{
    class_<HybridBox> hybrid_box_class("HybridBox",init<HybridBox>());
    hybrid_box_class.def(init<DiscreteLocation,Box>());
    hybrid_box_class.def("location",&HybridBox::location,return_value_policy<copy_const_reference>());
    hybrid_box_class.def("continuous_state_set",&HybridBox::continuous_state_set,return_value_policy<copy_const_reference>());
    hybrid_box_class.def(self_ns::str(self));
}

void export_hybrid_taylor_set()
{
    typedef HybridBasicSet<TaylorImageSet> HybridTaylorImageSet;

    class_<HybridTaylorImageSet> hybrid_taylor_set_class("HybridTaylorImageSet",init<HybridTaylorImageSet>());
    hybrid_taylor_set_class.def(init<DiscreteLocation,Box>());
    hybrid_taylor_set_class.def(init<DiscreteLocation,TaylorImageSet>());
    hybrid_taylor_set_class.def(self_ns::str(self));

    implicitly_convertible<HybridBox,HybridTaylorImageSet>();
}

void export_hybrid_constrained_image_set()
{
    typedef HybridBasicSet<ConstrainedImageSet> HybridConstrainedImageSet;

    class_<HybridConstrainedImageSet>
        hybrid_constrained_image_set_class("HybridConstrainedImageSet",init<HybridConstrainedImageSet>());
    hybrid_constrained_image_set_class.def(init<DiscreteLocation,Box>());
    //hybrid_constrained_image_set_class.def(init<DiscreteLocation,ConstrainedImageSet>());
    hybrid_constrained_image_set_class.def(self_ns::str(self));

    implicitly_convertible<HybridBox,HybridConstrainedImageSet>();
}




void geometry_submodule() {
    export_set_interface();
    export_point();
    export_box();
    export_zonotope();
    export_polytope();
    export_taylor_set();
    export_curve();

    export_affine_set();
    export_constrained_image_set();

    export_hybrid_box();
    export_hybrid_taylor_set();
    export_hybrid_constrained_image_set();

    typedef HybridBasicSet<TaylorImageSet> HybridTaylorImageSet;
    to_python< ListSet<TaylorImageSet> >();
    to_python< ListSet<HybridTaylorImageSet> >();
}

