/***************************************************************************
 *            dynamical_toolbox.cc
 *
 *  Copyright  2008  Pieter Collins
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
 
#include "macros.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "sparse_differential.h"
#include "function.h"
#include "box.h"
#include "dynamical_toolbox.h"

#include "approximate_taylor_model.h"



namespace Ariadne {

uint verbosity = 0;

class DegenerateCrossingException { };
class NonInvertibleFunctionException { };






/*
template<class Mdl>
typename DynamicalToolbox<Mdl>::ModelType
DynamicalToolbox<Mdl>::
model(const FunctionType& f) const
{
  Vector<I> domain=set.domain().position_vectors();
  Vector<R> midpoint=Ariadne::midpoint(domain);
  Vector<A> centre=Vector<A>(set.centre().position_vector());
  Matrix<A> const& generators=reinterpret_cast<Matrix<A>const&>(set.generators());
  return ApproximateTaylorModel<R>::affine(domain,midpoint,centre,generators,this->_spacial_order,this->_smoothness);
}
*/

  
/*
template<class Mdl>
typename DynamicalToolbox<Mdl>::ModelType
DynamicalToolbox<Mdl>::
set(const BoxType& box) const
{
  uint n=model.result_size();
  uint m=model.argument_size();
  Vector<I> set_centre(n);
  Matrix<R> set_generators(n,m);
  make_lpair(set_centre,set_generators) = affine_model(model);
  ES enclosure_set(Point<I>(set_centre),set_generators);
  return enclosure_set;
}
*/


  
template<class Mdl>
DynamicalToolbox<Mdl>::
DynamicalToolbox() 
  : _spacial_order(2),
    _temporal_order(4),
    _order(6),
    _smoothness(1)
{ 
}
  

template<class Mdl>
typename DynamicalToolbox<Mdl>::SetModelType
DynamicalToolbox<Mdl>::
integration_step(const FlowModelType& flow_model, 
                 const SetModelType& initial_set_model, 
                 const TimeType& integration_time) const
{
  uint ng=initial_set_model.argument_size();
  Mdl integration_time_model = Mdl::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),Vector<A>(1,A(integration_time)),_order,_smoothness);
  return this->integration_step(flow_model,initial_set_model,integration_time_model);                   }


template<class Mdl>
typename DynamicalToolbox<Mdl>::SetModelType
DynamicalToolbox<Mdl>::
integration_step(const FlowModelType& flow_model, 
                 const SetModelType& initial_set_model, 
                 const TimeModelType& integration_time_model) const
{
    Mdl set_step_model=join(initial_set_model, integration_time_model);
    ARIADNE_LOG(6,"set_step_model = "<<set_step_model<<"\n");
    Mdl final_set_model=compose(flow_model,set_step_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");

    return final_set_model;
}





template<class Mdl>
typename DynamicalToolbox<Mdl>::SetModelType
DynamicalToolbox<Mdl>::
reachability_step(const FlowModelType& flow_model, 
                  const SetModelType& initial_set_model, 
                  const TimeType& initial_time, 
                  const TimeType& final_time) const
{
  uint ng=initial_set_model.argument_size();
  Mdl initial_time_model = Mdl::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),
                                         Vector<A>(1,A(initial_time)),_order,_smoothness);
  Mdl final_time_model = Mdl::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),
                                       Vector<A>(1,A(final_time)),_order,_smoothness);
  return this->reachability_step(flow_model,initial_set_model,initial_time_model,final_time_model);
}


template<class Mdl>
typename DynamicalToolbox<Mdl>::ModelType
DynamicalToolbox<Mdl>::
reachability_step(const ModelType& flow_model, 
                  const ModelType& initial_set_model, 
                  const ModelType& initial_time_model, 
                  const ModelType& final_time_model) const
{
  // Compute the reachable set
  // Need an extra independent variable to represent time
  uint ng=initial_set_model.argument_size();
  
  ModelType expanded_initial_set_model=embed(initial_set_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
  ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
  ModelType expanded_initial_time_model=embed(initial_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
  ARIADNE_LOG(6,"expanded_initial_time_model="<<expanded_initial_time_model<<"\n");
  ModelType expanded_final_time_model=embed(final_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
  ARIADNE_LOG(6,"expanded_final_time_model="<<expanded_final_time_model<<"\n");
  
  ModelType time_interval_model=ModelType::affine(I(-1,1),R(0),A(0.5),A(0.5),_order,_smoothness);
  ARIADNE_LOG(6,"time_interval_time_model="<<time_interval_model<<"\n");
  ModelType expanded_time_interval_model=embed(time_interval_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),ng);
  ARIADNE_LOG(6,"expanded_time_interval_model="<<expanded_time_interval_model<<"\n");
  ModelType expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);
  ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
  ModelType expanded_timed_set_model=join(expanded_initial_set_model,expanded_reach_time_model);
  ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
  ModelType reach_set_model=compose(flow_model,expanded_timed_set_model);
  ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
  
  return reach_set_model;
}



template<class Mdl>
typename DynamicalToolbox<Mdl>::TimeModelType
DynamicalToolbox<Mdl>::
reachability_time(const TimeModelType& initial_time_model, 
                  const TimeModelType& final_time_model) const
{
  ARIADNE_ASSERT(initial_time_model.argument_size()==final_time_model.argument_size());
  uint ng=initial_time_model.argument_size();

  ModelType expanded_initial_time_model=embed(initial_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
  ARIADNE_LOG(6,"expanded_initial_time_model="<<expanded_initial_time_model<<"\n");
  ModelType expanded_final_time_model=embed(final_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
  ARIADNE_LOG(6,"expanded_final_time_model="<<expanded_final_time_model<<"\n");
  
  ModelType time_interval_model=ModelType::affine(I(-1,1),R(0),A(0.5),A(0.5),_order,_smoothness);
  ARIADNE_LOG(6,"time_interval_time_model="<<time_interval_model<<"\n");
  ModelType expanded_time_interval_model=embed(time_interval_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),ng);
  ARIADNE_LOG(6,"expanded_time_interval_model="<<expanded_time_interval_model<<"\n");
  ModelType expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);
  ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
  return expanded_reach_time_model;
}




// Compute the grazing time using bisections
template<class Mdl>
typename DynamicalToolbox<Mdl>::ModelType
DynamicalToolbox<Mdl>::
crossing_time(const ModelType& flow_model, 
              const ModelType& guard_model, 
              const ModelType& initial_set_model, 
              const RealType& minimum_time, 
              const RealType& maximum_time) const
{
  ARIADNE_ASSERT(minimum_time<=0);
  ARIADNE_ASSERT(maximum_time>=0);
  ModelType hitting_model=compose(guard_model,flow_model);
  ARIADNE_LOG(6,"hitting_model = "<<hitting_model<<"\n");
  ModelType free_hitting_time_model;
  try {
    free_hitting_time_model=implicit(hitting_model); 
  } catch(NonInvertibleFunctionException) {
    throw DegenerateCrossingException();
  }
  ARIADNE_LOG(6,"free_hitting_time_model = "<<free_hitting_time_model<<"\n");
  ModelType hitting_time_model=compose(free_hitting_time_model,initial_set_model);
  ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
  Interval hitting_time_range=hitting_time_model.range()[0];
  ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
  if(hitting_time_range.lower()<R(minimum_time) || hitting_time_range.upper()>R(maximum_time)) {
    throw DegenerateCrossingException();
  }

  return hitting_time_model;
}


// Compute the grazing time using bisections
template<class Mdl>
pair<Mdl,Mdl>
DynamicalToolbox<Mdl>::
touching_time_interval(const ModelType& flow_model, 
                       const ModelType& guard_model, 
                       const ModelType& initial_set_model, 
                       const RealType& minimum_time, 
                       const RealType& maximum_time) const
{
  ARIADNE_ASSERT(minimum_time<=0);
  ARIADNE_ASSERT(maximum_time>=0);
  ModelType final_set_model=this->integration_step(flow_model,initial_set_model,maximum_time);

  uint refinements=5;

  RealType lower_time=minimum_time;
  RealType upper_time=maximum_time;
  
  if(possibly(this->active(guard_model,initial_set_model))) {
    lower_time=minimum_time;
  } else {
    RealType min_lower_time=minimum_time;
    RealType max_lower_time=maximum_time;
    for(uint i=0; i!=refinements; ++i) {
      RealType new_lower_time=med_approx(min_lower_time,max_lower_time);
      if(possibly(this->active(guard_model,this->reachability_step(flow_model,initial_set_model,minimum_time,new_lower_time)))) {
        max_lower_time=new_lower_time;
      } else {
        min_lower_time=new_lower_time;
        lower_time=new_lower_time;
      }
    }
  }

  if(definitely(this->active(guard_model,initial_set_model)) || definitely(!this->active(guard_model,initial_set_model))) {
    upper_time=minimum_time;
  } else {
    RealType min_upper_time=minimum_time;
    RealType max_upper_time=maximum_time;
    for(uint i=0; i!=refinements; ++i) {
      RealType new_upper_time=med_approx(min_upper_time,max_upper_time);
      if(definitely(this->active(guard_model,this->integration_step(flow_model,initial_set_model,new_upper_time)))) {
        max_upper_time=new_upper_time;
        upper_time=new_upper_time;
      } else {
        min_upper_time=new_upper_time;
      }
    }
  }

  //return std::pair(R(lower_time),R(upper_time));
}






    


template<class Mdl>
tribool
DynamicalToolbox<Mdl>::
active(const ModelType& guard_model, const ModelType& set_model) const
{
  IntervalType range=compose(guard_model,set_model).range()[0];
  if(range.upper()<0) { 
    return false; 
  } else if(range.lower()>0) {
    return true;
  } else {
    return indeterminate;
  }
}

template<class Mdl>
tribool
DynamicalToolbox<Mdl>::
active(const FunctionType& guard_function, const ModelType& set_model) const
{
  ModelType guard_model(set_model.range(),guard_function,this->_spacial_order,this->_smoothness);
  return this->active(guard_model,set_model);
}



template class DynamicalToolbox<ApproximateTaylorModel>;

}  // namespace Ariadne
