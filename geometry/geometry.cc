/***************************************************************************
 *            geometry.cc
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

#include "numeric/numeric.h"
#include "config.h"

#include "geometry/geometry.h"
#include "utility/tuple.h"

namespace Ariadne {


tribool
separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps)
{
    ExactBox bb=make_exact_box(ls.bounding_box());
    if(bb.empty()) { return true; }
    return separated(ls,rs,bb,eps);
}


tribool
overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps)
{
    ExactBox bb=make_exact_box(ls.bounding_box());
    if(bb.empty()) { return false; }
    return overlap(ls,rs,bb,eps);
}


tribool
inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps)
{
    ExactBox bb=make_exact_box(ls.bounding_box());
    if(bb.empty()) { return true; }
    return inside(ls,rs,bb,eps);
}


tribool
overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBox& bx, const Float& eps)
{
    if(ls.separated(bx)) {
        return false;
    }
    if(rs.separated(bx)) {
        return false;
    }
    else if(rs.covers(bx)) {
        return true;
    }
    else if(bx.radius().raw()<eps) {
        return indeterminate;
    } else {
        ExactBox bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(ls.separated(bx1)) {
            return overlap(ls,rs,bx2,eps);
        } else if(ls.separated(bx2)) {
            return overlap(ls,rs,bx1,eps);
        } else {
            return overlap(ls,rs,bx1,eps) || overlap(ls,rs,bx2,eps);
        }
    }
}


tribool
inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBox& bx, const Float& eps)
{
    if(ls.separated(bx) || rs.separated(bx)) {
        return true;
    } else if(bx.radius().raw()<eps) {
        return indeterminate;
    } else {
        ExactBox bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(ls.separated(bx1)) {
            return inside(ls,rs,bx2,eps);
        } else if(ls.separated(bx2)) {
            return inside(ls,rs,bx1,eps);
        } else {
            return inside(ls,rs,bx1,eps) && inside(ls,rs,bx2,eps);
        }
    }
}


tribool
separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBox& bx, const Float& eps)
{
    if(ls.separated(bx) || rs.separated(bx)) {
        return true;
    } else if(bx.radius().raw()<eps) {
        return indeterminate;
    } else {
        ExactBox bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(ls.separated(bx1)) {
            return separated(ls,rs,bx2,eps);
        } else if(ls.separated(bx2)) {
            return separated(ls,rs,bx1,eps);
        } else {
            return separated(ls,rs,bx1,eps) && separated(ls,rs,bx2,eps);
        }
    }
}




tribool
overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const ExactBox& bx, const Float& eps)
{
    if(ovs.overlaps(bx)) {
        if(ops.covers(bx)) {
            return true;
        } else if(bx.radius().raw()<eps) {
            return indeterminate;
        } else {
            ExactBox bx1,bx2;
            make_lpair(bx1,bx2)=split(bx);
            if(overlap(ovs,ops,bx1,eps)) {
                return true;
            } else {
                return overlap(ovs,ops,bx2,eps);
            }
        }
    } else {
        return indeterminate;
    }
}


tribool
inside(const ClosedSetInterface& cls, const OpenSetInterface& ops, const ExactBox& bx, const Float& eps)
{
    if(cls.separated(bx) || ops.covers(bx)) {
        return true;
    } else if(bx.radius().raw()<eps) {
        return indeterminate;
    } else {
        ExactBox bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(inside(cls,ops,bx1,eps)) {
            return inside(cls,ops,bx2,eps);
        } else {
            return indeterminate;
        }
    }
}


tribool
separated(const ClosedSetInterface& cls1, const ClosedSetInterface& cls2, const ExactBox& bx, const Float& eps)
{
    if(cls1.separated(bx) || cls2.separated(bx)) {
        return true;
    } else if(bx.radius().raw()<eps) {
        return indeterminate;
    } else {
        ExactBox bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(separated(cls1,cls2,bx1,eps)) {
            return separated(cls1,cls2,bx2,eps);
        } else {
            return indeterminate;
        }
    }
}




} // namespace Ariadne
