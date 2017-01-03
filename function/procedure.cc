/***************************************************************************
 *            procedure.cc
 *
 *  Copyright 2016  Pieter Collins
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

#include "procedure.h"
#include "procedure.tcc"

#include "algebra/differential.h"
#include "algebra/graded.h"

namespace Ariadne {

template class Procedure<ApproximateNumber>;
template class Procedure<ValidatedNumber>;

template class Vector<Procedure<ApproximateNumber>>;
template class Vector<Procedure<ValidatedNumber>>;

template Void _execute(List<Float64Bounds>& v, const List<ProcedureInstruction>& p, const List<ValidatedNumber>& c, const Vector<Float64Bounds>& x);
template Void _execute(List<Graded<Differential<Float64Bounds>>>& v, const List<ProcedureInstruction>& p, const List<ValidatedNumber>& c, const Vector<Graded<Differential<Float64Bounds>>>& x);

inline
Void restrict(UpperIntervalType& r, const UpperIntervalType& x) {
    r.set_lower(max(r.lower(),x.lower()));
    r.set_upper(min(r.upper(),x.upper()));
};

Void simple_hull_reduce(UpperBoxType& dom, const ValidatedProcedure& f, ExactIntervalType codom)
{
    const List<ProcedureInstruction>& p=f._instructions;
    const List<ValidatedNumber>& c=f._constants;
    List<UpperIntervalType> t(p.size());
    _execute(t,p,c,dom);
    restrict(t.back(),codom);
    _backpropagate(dom,t,p,c);
}

Void simple_hull_reduce(UpperBoxType& dom, const Vector<ValidatedProcedure>& f, ExactBoxType codom)
{
    const List<ProcedureInstruction>& p=f._instructions;
    const List<ValidatedNumber>& c=f._constants;
    List<UpperIntervalType> t(p.size());

    ARIADNE_ASSERT(codom.size()==f._results.size());

    _execute(t,p,c,dom);
    for(SizeType i=0; i!=codom.size(); ++i) {
        restrict(t[f._results[i]],codom[i]);
    }
    _backpropagate(dom,t,p,c);
}

} // namespace Ariadne