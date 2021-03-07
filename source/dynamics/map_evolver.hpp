/***************************************************************************
 *            dynamics/map_evolver.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file dynamics/map_evolver.hpp
 *  \brief Evolver for map systems.
 */

#ifndef ARIADNE_MAP_EVOLVER_HPP
#define ARIADNE_MAP_EVOLVER_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>


#include "utility/tuple.hpp"

#include "dynamics/map.hpp"
#include "function/function_interface.hpp"
#include "configuration/configuration_interface.hpp"
#include "dynamics/evolver_base.hpp"

#include "output/logging.hpp"

namespace Ariadne {

template<class Sys, class BS, class TRM> class Evolver;

class Enclosure;
class IteratedMap;
template<class ES> class Orbit;

class MapEvolverConfiguration;
class EvolutionProfiler;


/*! \brief A class for computing the evolution of an iterated map.
 */
class MapEvolver
    : public EvolverBase< IteratedMap, LabelledEnclosure, Integer>
{
  public:
    typedef MapEvolverConfiguration ConfigurationType;
    typedef IteratedMap SystemType;
    typedef Integer TimeType;
    typedef Integer TerminationType;
    typedef LabelledEnclosure EnclosureType;
    typedef Pair<TerminationType, EnclosureType> TimedEnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef ValidatedFunctionModelDPFactory FunctionFactoryType;
  public:

    //! \brief Default constructor.
    MapEvolver(const SystemType& system);

    /*! \brief Make a dynamically-allocated copy. */
    MapEvolver* clone() const { return new MapEvolver(*this); }

    /* \brief Get the internal system. */
    virtual const SystemType& system() const { return *_sys_ptr; }

    //! \brief Make an enclosure from a user set.
    EnclosureType enclosure(RealBox const&) const;
    EnclosureType enclosure(RealBox const&, EnclosureConfiguration const&) const;
    EnclosureType enclosure(RealVariablesBox const&) const;
    EnclosureType enclosure(RealVariablesBox const&, EnclosureConfiguration const&) const;
    //! \brief Make an enclosure from a computed box set.
    EnclosureType enclosure(ExactBoxType const&) const;
    EnclosureType enclosure(ExactBoxType const&, EnclosureConfiguration const&) const;

    //!@{
    //! \name Configuration for the class.

    //! \brief A reference to the configuration.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //! \brief The class which constructs functions for the enclosures.
    const FunctionFactoryType function_factory() const;

    //!@}


    //!@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;


    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,false);
        return final; }

    //! \brief Compute an approximation to the evolution set under upper semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,true);
        return intermediate; }


    Pair<EnclosureListType,EnclosureListType> reach_evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,true); return std::make_pair(reachable,final); }


  protected:
    virtual Void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const EnclosureType& initial, const TerminationType& termination,
                            Semantics semantics, Bool reach) const;

    virtual Void _evolution_step(List< TimedEnclosureType >& working_sets,
                                 EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                 const TimedEnclosureType& current_set, const TerminationType& termination,
                                 Semantics semantics, Bool reach) const;

  private:
    std::shared_ptr< SystemType > _sys_ptr;
    //std::shared_ptr< EvolutionProfiler >  _profiler;
    std::shared_ptr< ConfigurationType > _configuration;
};


//! \brief Configuration for a MapEvolver, essentially to control accuracy of evolution.
class MapEvolverConfiguration : public ConfigurationInterface
{
  public:
    typedef double RealType;

    //! \brief Default constructor gives reasonable values.
    MapEvolverConfiguration();

    virtual ~MapEvolverConfiguration() = default;

  private:

    //! \brief The maximum allowable radius of a basic set.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_enclosure_radius;

    //! \brief Allow subdivisions in upper evolution.
    Bool _enable_subdivisions;

    //! \brief Allow premature termination of lower evolution.
    Bool _enable_premature_termination;

  public:

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    Void set_maximum_enclosure_radius(const RealType value) { _maximum_enclosure_radius = value; }

    const Bool& enable_subdivisions() const { return _enable_subdivisions; }
    Void set_enable_subdivisions(const Bool value) { _enable_subdivisions = value; }

    const Bool& enable_premature_termination() const { return _enable_premature_termination; }
    Void set_enable_premature_termination(const Bool value) { _enable_premature_termination = value; }

  public:

    virtual OutputStream& _write(OutputStream& os) const;

};


} // namespace Ariadne

#endif // ARIADNE_MAP_EVOLVER_HPP
