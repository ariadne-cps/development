/***************************************************************************
 *            hybrid_set.h
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

/*! \file hybrid_set.h
 *  \brief Sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_H
#define ARIADNE_HYBRID_SET_H

#include <map>

#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include <memory>

#include "macros.h"
#include "stlio.h"
#include "container.h"
#include "function_set.h"
#include "list_set.h"
#include "grid_set.h"
#include "curve.h"

#include "expression_set.h"

#include "hybrid_set_interface.h"
#include "hybrid_space.h"
#include "hybrid_grid.h"
#include "point.h"
#include "box.h"

#ifdef ARIADNE_ENABLE_SERIALIZATION
#include "serialization.h"
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#include "graphics_interface.h"

namespace Ariadne {

class HybridGridTreeSet;

class RealVariableInterval;
class RealVariablesBox;

typedef Map< DiscreteLocation, Vector<Float> > HybridFloatVector;

template<class HBS> class HybridBasicSetExpression { };
template<class HDS> class HybridDenotableSetExpression { };


template<class SET> struct is_basic_set { };

// Declare template specialisation for hybrid list set
template<class ES> class HybridBasicSet;
template<class ES> class ListSet< HybridBasicSet<ES> >;

template<class DS, class HBS> class HybridSetConstIterator;

class HybridRealExpressionBoundedConstraintSet;
typedef HybridRealExpressionBoundedConstraintSet HybridExpressionSet;

class HybridRealBoundedConstraintSet;
typedef HybridRealBoundedConstraintSet HybridSet;

class HybridRealExpressionBoundedConstraintSet
    : public HybridDrawableInterface
{
    DiscreteLocation _location;
    Map< RealVariable, RealIntervalSet> _bounds;
    List< ContinuousPredicate > _constraints;
  public:
    HybridRealExpressionBoundedConstraintSet() { }
    HybridRealExpressionBoundedConstraintSet(const DiscreteLocation& loc, const List<RealVariableInterval>& bnd, const List<ContinuousPredicate>& cnstr = List<ContinuousPredicate>());
    HybridRealExpressionBoundedConstraintSet(const DiscreteLocation& loc, const RealExpressionBoundedConstraintSet& set);
    HybridRealExpressionBoundedConstraintSet* clone() const { return new HybridRealExpressionBoundedConstraintSet(*this); }
    DiscreteLocation const& location() const { return this->_location; }
    Set<RealVariable> variables() const { return this->_bounds.keys(); };
    Map<RealVariable,RealIntervalSet> const& bounds() const { return this->_bounds; };
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; };
    RealBoundedConstraintSet continuous_set(const RealSpace&) const;
    void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const;
};
OutputStream& operator<<(OutputStream& os, const HybridRealExpressionBoundedConstraintSet& hs);

typedef HybridRealExpressionBoundedConstraintSet HybridExpressionSet;

class HybridPoint
    : public Tuple<DiscreteLocation,RealSpace,Point>
{
  public:
    HybridPoint() : Tuple<DiscreteLocation,RealSpace,Point>(DiscreteLocation(),RealSpace(),Point()) { }
    HybridPoint(const DiscreteLocation& q, const RealSpace& spc, const Point& pt) : Tuple<DiscreteLocation,RealSpace,Point>(q,spc,pt) { }
    HybridPoint(const DiscreteLocation& q, const Map<Identifier,Float>& val);
    HybridPoint(const DiscreteLocation& q, const Map<Identifier,Real>& val);
    HybridPoint(const DiscreteLocation& q, const List< Assignment<RealVariable,Real> >& val);
    const DiscreteLocation& location() const { return this->first; }
    const RealSpace& space() const { return this->second; }
    const Point& point() const { return this->third; }
    const Point& real_values() const { return this->third; }
    Map<RealVariable,Float> values() const;
};

template<class BS>
class HybridBasicSet
    : public Tuple<DiscreteLocation,RealSpace,BS>
{
  public:
    typedef BS ContinuousSetType;
    HybridBasicSet() : Tuple<DiscreteLocation,RealSpace,BS>(DiscreteLocation(),RealSpace(),BS()) { }
    HybridBasicSet(const DiscreteLocation& q, const RealSpace& spc, const BS& bs) : Tuple<DiscreteLocation,RealSpace,BS>(q,spc,bs) { }
    const DiscreteLocation& location() const { return this->first; }
    const RealSpace& space() const { return this->second; }
    const ContinuousSetType& continuous_set() const { return this->third; }
};

typedef HybridBasicSet<Box> HybridBox;

template<> class HybridBasicSet<Box>
    : public HybridDrawableInterface
{
    DiscreteLocation _location;
    RealSpace _space;
    Box _set;
  public:
    HybridBasicSet() : _location(), _space(), _set() { }
    HybridBasicSet(const DiscreteLocation& q, const RealSpace& spc, const Box& bx) : _location(q), _space(spc), _set(bx) { }
    HybridBasicSet(const DiscreteLocation& q, const List<RealVariableInterval>& bnds);
    const DiscreteLocation& location() const { return this->_location; }
    const RealSpace& space() const { return this->_space; }
    const Box& continuous_set() const { return this->_set; }
    const Interval& operator[](const RealVariable& v) const { return this->_set[this->_space.index(v)]; }
    HybridBasicSet<Box>* clone() const { return new HybridBasicSet<Box>(*this); }
    void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const;
    friend std::ostream& operator<<(std::ostream& os, const HybridBox& hbx) {
        return os << "(" << hbx.location() << ", " << hbx.space() << ", " << hbx.continuous_set() << ")"; }
};

template<class ES> inline
bool operator==(const HybridBasicSet<ES> hset1, const HybridBasicSet<ES>& hset2) {
    return hset1.location()==hset2.location() && hset1.continuous_set() == hset2.continuous_set();
}


typedef List<RealVariable> RealVariables;

class HybridRealBoundedConstraintSet
    : public virtual HybridSetInterface
    , public virtual HybridDrawableInterface
{
    Map<DiscreteLocation, RealExpressionBoundedConstraintSet> _sets;
    mutable Map<DiscreteLocation,Pair<RealSpace,RealBoundedConstraintSet>> _euclidean_sets;
  public:
    HybridRealBoundedConstraintSet();
    HybridRealBoundedConstraintSet(const HybridBox& bx);
    HybridRealBoundedConstraintSet(const DiscreteLocation& loc,
                                   const RealVariablesBox& bx);
    HybridRealBoundedConstraintSet(const DiscreteLocation& loc,
                                   const List<RealVariableInterval>& bnd,
                                   const List<ContinuousPredicate>& cnstr = List<ContinuousPredicate>());

    DiscreteLocation location() const;

    virtual HybridRealBoundedConstraintSet* clone() const override;
    //virtual HybridSpace space() const;
    virtual tribool overlaps(const HybridBox& bx) const override;
    virtual tribool inside(const HybridBoxes& bx) const override;

    virtual tribool separated(const HybridBox& bx) const override;
    virtual tribool covers(const HybridBox& bx) const override;
    virtual RealBoundedConstraintSet const& euclidean_set(DiscreteLocation loc, RealSpace spc) const override;
    virtual Set<DiscreteLocation> locations() const override;
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    virtual HybridBoxes bounding_box() const override;
    virtual std::ostream& write(std::ostream& os) const override;
    virtual void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
};


template<class ES> class HybridListSetConstIterator;

template<class ES>
class HybridListSetConstIterator
    : public boost::iterator_facade<
                HybridListSetConstIterator<ES>,
                HybridBasicSet<ES>,
                boost::forward_traversal_tag,
                HybridBasicSet<ES> const&
             >

{
  public:
    typedef HybridBasicSet<ES> const& reference;
  public:
    HybridListSetConstIterator(const Map<DiscreteLocation,Pair<RealSpace,ListSet<ES> > >& map, bool);
    bool equal(const HybridListSetConstIterator<ES>&) const;
    const HybridBasicSet<ES>& dereference() const;
    void increment();
  private:
    void _increment_loc();
  private:
    typedef typename Map< DiscreteLocation, Pair<RealSpace,ListSet<ES> > >::const_iterator LocationsIterator;
    typedef typename ListSet<ES>::const_iterator BasicSetIterator;
    LocationsIterator _loc_begin;
    LocationsIterator _loc_end;
    LocationsIterator _loc_iter;
    BasicSetIterator _bs_iter;
    mutable HybridBasicSet<ES> hybrid_set;
};


template<class ES> inline
HybridListSetConstIterator<ES>::
HybridListSetConstIterator(const Map<DiscreteLocation,Pair<RealSpace,ListSet<ES> > >& map, bool end)
    : _loc_begin(map.begin()),
      _loc_end(map.end()),
      _loc_iter(end?_loc_end:_loc_begin)
{
    if(_loc_iter!=_loc_end) {
        _bs_iter=_loc_iter->second.second.begin();
        this->_increment_loc();
    }
}


template<class ES> inline
bool
HybridListSetConstIterator<ES>::equal(const HybridListSetConstIterator<ES>& other) const
{
    return this->_loc_iter==other._loc_iter && (this->_loc_iter==this->_loc_end || this->_bs_iter==other._bs_iter);
}


template<class ES> inline
HybridBasicSet<ES> const&
HybridListSetConstIterator<ES>::dereference() const
{
    this->hybrid_set=HybridBasicSet<ES>(_loc_iter->first,_loc_iter->second.first,*this->_bs_iter);
    return this->hybrid_set;
}


template<class ES> inline
void
HybridListSetConstIterator<ES>::increment()
{
    ++this->_bs_iter;
    this->_increment_loc();
}

template<class ES> inline
void
HybridListSetConstIterator<ES>::_increment_loc()
{
    while(_bs_iter==_loc_iter->second.second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.second.begin();
    }
}


template<class ES> class HybridListSet;
template<class ES> std::ostream& operator<<(std::ostream& os, const HybridListSet<ES>& hls);

//! \ingroup HybridModule
//! A set comprising a %ListSet in each location.
template<class ES>
class HybridListSet
{
    Map<DiscreteLocation, Pair< RealSpace, ListSet<ES> > > _locations;
  public:
    typedef typename Map<DiscreteLocation,Pair< RealSpace, ListSet<ES> > >::const_iterator locations_const_iterator;
    typedef typename Map<DiscreteLocation,Pair< RealSpace, ListSet<ES> > >::iterator locations_iterator;
    typedef HybridListSetConstIterator<ES> const_iterator;

    HybridListSet() { }
    HybridListSet(const HybridBasicSet<ES>& hes) { this->adjoin(hes); }

    const_iterator begin() const {
        return const_iterator(this->_locations,false); }
    const_iterator end() const {
        return const_iterator(this->_locations,true); }

    locations_const_iterator locations_begin() const {
        return this->_locations.begin(); }
    locations_const_iterator locations_end() const {
        return this->_locations.end(); }

    //! \brief Returns the number of basic hybrid sets forming this object.
    size_t size() const {
        size_t s = 0;
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            s += _loc_iter->second.second.size();
        }
        return s;
    }

    //using std::map< DiscreteLocation,ListSet<ES> >::insert;

    //using std::map<DiscreteLocation,ListSet<ES> >::operator[];
    //ListSet<ES>& operator[](const DiscreteLocation& q) {
    //    return this->_locations[q].second; }
    const ListSet<ES>& operator[](const DiscreteLocation& q) const {
        ARIADNE_ASSERT_MSG(this->find(q)!=this->locations_end(),(*this)<<" has no location "<<q);
        return this->find(q)->second->second; }

    void insert(const DiscreteLocation& loc, const RealSpace& spc) {
        this->_locations.insert(loc, make_pair(spc,ListSet<ES>())); }
    void adjoin(const DiscreteLocation& loc, const RealSpace& spc, const ES& es) {
        locations_iterator loc_iter=this->_locations.find(loc);
        if(loc_iter!=this->_locations.end()) {
            ARIADNE_ASSERT_MSG(loc_iter->second.first==spc,
                               "Space "<<loc_iter->second.first<<" of location "<<loc_iter->first<<" differs from "<<spc); }
        else { this->insert(loc,spc); loc_iter=this->_locations.find(loc); }
        loc_iter->second.second.adjoin(es); }
    void adjoin(const DiscreteLocation& loc, const ES& es) {
        ARIADNE_ASSERT_MSG(this->_locations.find(loc)!=this->_locations.end(),(*this)<<" has no location "<<loc);
        this->_locations[loc].second.adjoin(es); }
    void adjoin(const HybridBasicSet<ES>& hes) {
        this->adjoin(hes.location(),hes.space(),hes.continuous_set()); }
    void adjoin(const HybridListSet<ES>& hls) {
        for(locations_const_iterator _loc_iter=hls.locations_begin();
            _loc_iter!=hls.locations_end(); ++_loc_iter) {
            (*this)[_loc_iter->first].adjoin(_loc_iter->second); } }

    HybridListSet<Box> bounding_boxes() const {
        HybridListSet<Box> result;
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            result[_loc_iter->first]=_loc_iter->second.bounding_boxes(); }
        return result; }

    friend
    std::ostream& operator<<(std::ostream& os, const HybridListSet<ES>& hls) {
        return os << hls._locations; }
};

template<class ES>
class ListSet< HybridBasicSet<ES> >
    : public HybridListSet<ES>
{
  public:
    ListSet() { }
    ListSet(const HybridBasicSet<ES>& hes) { this->adjoin(hes); }
};


template<class ES>
std::ostream&
operator<<(std::ostream& os,
           const ListSet< HybridBasicSet<ES> >& ls) {
    const std::map< DiscreteLocation, ListSet<ES> >& hls=ls;
    return os << "HybridListSet" << hls;
}




class HybridGridCell
    : public HybridBasicSet<GridCell>
{
  public:
    HybridGridCell()
        : HybridBasicSet<GridCell>() { }
    HybridGridCell(DiscreteLocation q,const RealSpace& s,const GridCell& gc)
        : HybridBasicSet<GridCell>(q,s,gc) { }
    HybridBox box() const { return HybridBox(this->first,this->second,this->third.box()); }
};

class HybridGridTreeSet;



template<class HDS1, class HDS2>
void adjoin_denotable_set(HDS1& hds1, const HDS2 hds2) {
    for(typename HDS2::locations_const_iterator loc2_iter=
            hds2.locations_begin(); loc2_iter!=hds2.locations_end(); ++loc2_iter)
        {
            typename HDS1::locations_iterator loc1_iter=hds1.find(loc2_iter->first);
            if(loc1_iter==hds1.locations_end()) {
                hds1.insert(make_pair(loc2_iter->first,
                                      typename HDS2::value_type::second_type(loc2_iter->second)));
            } else {
                loc1_iter->second.adjoin(loc2_iter->second);
            }
        }
}


template<class DS, class HBS>
class HybridSetConstIterator
    : public boost::iterator_facade<HybridSetConstIterator<DS,HBS>,
                                    HBS,
                                    boost::forward_traversal_tag,
                                    HBS const&
                                    >

{
  public:
    typedef HBS const& reference;
  public:
    HybridSetConstIterator(const std::map<DiscreteLocation,DS>&, const HybridSpace& hspc, bool);
    bool equal(const HybridSetConstIterator<DS,HBS>&) const;
    const HBS& dereference() const;
    void increment();
  private:
    void increment_loc();
  private:
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_begin;
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_end;
    typename std::map< DiscreteLocation,DS>::const_iterator _loc_iter;
    typename DS::const_iterator _bs_iter;
    HybridSpace hspc;
    mutable HBS hybrid_set;
};



//! \ingroup HybridModule
//! A set comprising a %GridTreeSet in each location.
class HybridGridTreeSet
    : public HybridDrawableInterface
{
  public:
    HybridGrid _hgrid;
    Map<DiscreteLocation,GridTreeSet> _map;
  public:
    typedef Map<DiscreteLocation,GridTreeSet>::iterator locations_iterator;
    typedef Map<DiscreteLocation,GridTreeSet>::const_iterator locations_const_iterator;
    typedef HybridSetConstIterator<GridTreeSet,HybridGridCell> const_iterator;
  public:
    //!
    locations_iterator locations_begin() { return this->_map.begin(); }
    //!
    locations_iterator locations_end() { return this->_map.end(); }
    //!
    locations_const_iterator locations_begin() const { return this->_map.begin(); }
    //!
    locations_const_iterator locations_end() const { return this->_map.end(); }
    //!
    const_iterator begin() const { return const_iterator(this->_map,this->_hgrid.space(),false); }
    //!
    const_iterator end() const { return const_iterator(this->_map,this->_hgrid.space(),true); }
  public:
    //! Construct from a grid.
    HybridGridTreeSet(const HybridGrid& hgrid) : _hgrid(hgrid), _map() { }

    //!
    HybridGrid grid() const { return this->_hgrid; }

    //! Test if \a q is a location of the set i.e. corresponds to a valid location of the underlying hybrid space.
    bool has_location(DiscreteLocation q) const { return _hgrid.has_location(q); }
    //! Test if \a q is a nontrivial location of the set i.e. contained in the map of <DiscreteLocation,GridTreeSet> pairs.
    bool nontrivial_location(DiscreteLocation q) const { return _map.has_key(q); }
    //! The continuous state space corresponding to location \a q.
    RealSpace space(DiscreteLocation q) const { return _hgrid.space(q); }
    //! The continuous state space corresponding to location \a q.
    const GridTreeSet& continuous_set(DiscreteLocation q) const { return _map[q]; }
    //! The continuous state space corresponding to location \a q.
    const GridTreeSet& euclidean_set(DiscreteLocation q, const RealSpace& s) const {
        ARIADNE_ASSERT_MSG(s==this->space(q),"Variable ordering in HybridGridTreeSet location "<<q<<" is "<<this->space(q)<<", "
                                             "which does not match requested ordering "<<q);
        return _map[q]; }

    //!
    void insert(DiscreteLocation q, const GridTreeSet& gts) {
        this->_map.insert(q,gts); }
    void insert(std::pair<DiscreteLocation,GridTreeSet>& qgts) {
        this->_map.insert(qgts); }

    //!
    void adjoin(DiscreteLocation q, const GridCell& c) {
        this->_provide_location(q).adjoin(c); }

    //!
    void adjoin(const HybridGridCell& hgc) {
        this->_provide_location(hgc.first).adjoin(hgc.third); }

    //!
    void adjoin(const ListSet<HybridGridCell>& hgcls) {
        for(ListSet<HybridGridCell>::const_iterator iter=hgcls.begin(); iter!=hgcls.end(); ++iter) {
            this ->adjoin(*iter); } }

    //!
    void adjoin(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).adjoin(_loc_iter->second); } }

    //!
    void remove(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).remove(_loc_iter->second); } }

    //!
    void restrict(const HybridGridTreeSet& hgts) {
        for(HybridGridTreeSet::locations_const_iterator _loc_iter=hgts.locations_begin(); _loc_iter!=hgts.locations_end(); ++_loc_iter) {
            this->_provide_location(_loc_iter->first).restrict(_loc_iter->second); } }

    //!
    void restrict_to_height(uint h) {
        for(locations_iterator _loc_iter=locations_begin(); _loc_iter!=locations_end(); ++_loc_iter) {
            _loc_iter->second.restrict_to_height(h); } }

    //!
    void adjoin_inner_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator _loc_iter=hbxs.begin();
                _loc_iter!=hbxs.end(); ++_loc_iter) {
            DiscreteLocation loc=_loc_iter->first;
            RealSpace spc=this->space(loc);
            Box const& bx = _loc_iter->second.euclidean_set(spc.variables());
            this->_provide_location(loc).adjoin_inner_approximation(bx,bx,depth); } }

    //!
    void adjoin_lower_approximation(const HybridOvertSetInterface& hs, const int height, const int depth) {
        Set<DiscreteLocation> hlocs=dynamic_cast<const HybridBoundedSetInterface&>(hs).locations();
        for(Set<DiscreteLocation>::const_iterator _loc_iter=hlocs.begin();
                _loc_iter!=hlocs.end(); ++_loc_iter) {
            DiscreteLocation loc=*_loc_iter;
            RealSpace spc=this->space(loc);
            this->_provide_location(loc).adjoin_lower_approximation(hs.euclidean_set(loc,spc),height,depth); } }

    //!
    void adjoin_outer_approximation(const HybridCompactSetInterface& hs, const int depth) {
        Set<DiscreteLocation> hlocs=hs.locations();
        for(Set<DiscreteLocation>::const_iterator _loc_iter=hlocs.begin();
                _loc_iter!=hlocs.end(); ++_loc_iter) {
            DiscreteLocation loc=*_loc_iter;
            RealSpace spc=this->space(loc);
            this->_provide_location(loc).adjoin_outer_approximation(hs.euclidean_set(loc,spc),depth); } }

    //!
    void adjoin_outer_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator _loc_iter=hbxs.begin();
                _loc_iter!=hbxs.end(); ++_loc_iter) {
            DiscreteLocation loc=_loc_iter->first;
            RealSpace spc=this->space(loc);
            Box const& bx = _loc_iter->second.euclidean_set(spc.variables());
            this->_provide_location(_loc_iter->first).adjoin_outer_approximation(bx,depth); } }

    //!
    template<class S> void adjoin_outer_approximation(DiscreteLocation q, const S& s) {
        this->_provide_location(q).adjoin_outer_approximation(s); }

    //!
    GridTreeSet& operator[](DiscreteLocation q) {
        return this->_provide_location(q);
    }

    //!
    const GridTreeSet& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT_MSG(this->has_location(q),"q="<<q);
        return const_cast<HybridGridTreeSet*>(this)->_provide_location(q);
    }

    //!
    bool empty() const {
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            if(!_loc_iter->second.empty()) { return false; } }
        return true; }

    //!
    size_t size() const {
        size_t result=0;
        for(locations_const_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            result+=_loc_iter->second.size(); }
        return result; }

    //!
    HybridListSet<Box> boxes() const {
        HybridListSet<Box> result;
        for(const_iterator iter=this->begin();
            iter!=this->end(); ++iter) {
            result.adjoin(iter->first,iter->third.box()); }
        return result; }

    //!
    void mince(int depth) {
        for(locations_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            _loc_iter->second.mince(depth); } }

    //!
    void recombine() {
        for(locations_iterator _loc_iter=this->locations_begin();
            _loc_iter!=this->locations_end(); ++_loc_iter) {
            _loc_iter->second.recombine(); } }
  public:
    //@{ \name HybridSetInterface methods

    //!
    HybridGridTreeSet* clone() const { return new HybridGridTreeSet(*this); }

    //!
    HybridSpace space() const { return this->grid().space(); }

    //!
    tribool separated(const HybridBox& hbx) const {
        locations_const_iterator _loc_iter = this->_map.find( hbx.location() );
        return _loc_iter != this->locations_end() || _loc_iter->second.separated( hbx.continuous_set() );
    }

    //!
    tribool overlaps(const HybridBox& hbx) const {
        locations_const_iterator _loc_iter = this->_map.find( hbx.location() );
        return _loc_iter != this->locations_end() && _loc_iter->second.overlaps( hbx.continuous_set() );
    }

    //!
    tribool covers(const HybridBox& hbx) const {
        locations_const_iterator _loc_iter=this->_map.find(hbx.location());
        return _loc_iter!=this->locations_end() && _loc_iter->second.covers( hbx.continuous_set() );
    }

    //!
    tribool inside(const HybridBoxes& hbx) const  {
        for( locations_const_iterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
            if( !_loc_iter->second.empty() ) {
                DiscreteLocation const& loc = _loc_iter->first;
                RealSpace spc=this->space(loc);
                return this->continuous_set(loc).inside(hbx.euclidean_set(loc,spc));
            }
        }
        return true;
    }

    //!
    HybridBoxes bounding_box() const {
        HybridBoxes result;
        for( locations_const_iterator _loc_iter = this->locations_begin(); _loc_iter != this->locations_end(); ++_loc_iter ) {
            if( !_loc_iter->second.empty() ) {
                DiscreteLocation const& loc = _loc_iter->first;
                RealSpace const& spc=this->space(loc);
                result.insert(loc,spc,_loc_iter->second.bounding_box());
            }
        }
        return result;
    }

    //!
    std::ostream& write(std::ostream& os) const {
        return os << this->_map;
    }

    //!
    void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d&v) const;

    //!
    friend inline std::ostream&
    operator<<(std::ostream& os, const HybridGridTreeSet& hgts) {
        return os << "HybridGridTreeSet(" << hgts._map << ")"; }

  private:
    GridTreeSet& _provide_location(const DiscreteLocation& q) {
        std::map<DiscreteLocation,GridTreeSet>::iterator iter=this->_map.find(q);
        if(iter==this->_map.end()) {
            this->_map.insert(std::make_pair(q,GridTreeSet(this->_hgrid[q])));
            iter=this->_map.find(q);
        }
        return iter->second; }
};




template<class DS, class HBS> inline
HybridSetConstIterator<DS,HBS>::
HybridSetConstIterator(const std::map<DiscreteLocation,DS>& map, const HybridSpace& spc, bool end)
    : _loc_begin(map.begin()),
      _loc_end(map.end()),
      _loc_iter(end?_loc_end:_loc_begin),
      hspc(spc)
{
    if(_loc_iter!=_loc_end) {
        _bs_iter=_loc_iter->second.begin();
        this->increment_loc();
    }
}


template<class DS, class HBS> inline
bool
HybridSetConstIterator<DS,HBS>::equal(const HybridSetConstIterator<DS,HBS>& other) const
{
    return this->_loc_iter==other._loc_iter && (this->_loc_iter==this->_loc_end || this->_bs_iter==other._bs_iter);
}


template<class DS, class HBS> inline
HBS const&
HybridSetConstIterator<DS,HBS>::dereference() const
{
    this->hybrid_set=HBS(_loc_iter->first,this->hspc[_loc_iter->first],*this->_bs_iter);
    return this->hybrid_set;
}


template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment()
{
    ++this->_bs_iter;
    this->increment_loc();
}

template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment_loc()
{
    while(_bs_iter==_loc_iter->second.end()) {
        ++_loc_iter;
        if(_loc_iter==_loc_end) { return; }
        _bs_iter=_loc_iter->second.begin();
    }
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> void serialize(A& archive, HybridGridTreeSet& set, const unsigned int version) {
      archive & static_cast<std::map<DiscreteLocation,GridTreeSet>&>(set); }
#endif /* ARIADNE_ENABLE_SERIALIZATION */


} // namespace Ariadne

#ifdef ARIADNE_ENABLE_SERIALIZATION
  namespace boost { namespace serialization {
  template<class A> void serialize(A& archive, const Ariadne::HybridGridTreeSet& set, const uint version);
  template<class A> void serialize(A& archive, const Ariadne::DiscreteLocation& state, const uint version);
  }}
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#endif // ARIADNE_HYBRID_SET_H
