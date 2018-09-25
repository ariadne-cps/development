/***************************************************************************
 *            numeric/paradigm.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/paradigm.hpp
 *  \brief
 */

#ifndef ARIADNE_PARADIGM_HPP
#define ARIADNE_PARADIGM_HPP

#include <cstdint>
#include "../utility/metaprogramming.hpp"

namespace Ariadne {

class ParadigmError { };

typedef std::uint16_t  ParadigmCodeType;

enum class ParadigmCode : ParadigmCodeType {
    APPROXIMATE=1,
    VALIDATED=2,
    EFFECTIVE=3,
    EXACT=4
};

//! \ingroup ParadigmSubModule
//! \brief The <em>computational paradigm</em> supported by the object.
//! User paradigms are ExactTag, EffectiveTag, ValidatedTag,ValidatedBoundedTag, ValidatedUpperTag, ValidatedLowerTag or ApproximateTag.
//! Internal paradigms are BuiltinTag and RawTag.
template<class T> using Paradigm = typename T::Paradigm;

//! \ingroup ParadigmSubModule
//! \brief The <em>computational paradigm</em> supported by the object. Equivalent to Paradigm<T>.
template<class T> using ParadigmTag = typename T::Paradigm;

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object is of a builtin type. Such objects should be converted to %Ariadne internal types before use.
struct BuiltinTag { };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object decribes raw data. Such objects should not be used in high-level code, as they probably do not provide safe guarantees on their values.
struct RawTag { };


//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object provides an approximation to a quantity with no guarantees on the error.
struct ApproximateTag { static constexpr ParadigmCode code() { return ParadigmCode::APPROXIMATE; } };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents an approximation to a quantity with a bound on the error in some metric, or lower and upper bounds on a quantity.
struct ValidatedTag : public ApproximateTag { static constexpr ParadigmCode code() { return ParadigmCode::VALIDATED; } };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, but equality is undecidable.
struct EffectiveTag : public ValidatedTag { static constexpr ParadigmCode code() { return ParadigmCode::EFFECTIVE; } };

//! \ingroup ParadigmSubModule
//! \brief A tag meaning that the object represents a quantity exactly, and equality is decidable.
//! Only available for discrete types i.e. elements of countable spaces, or for computational types such as floating-point numbers.
struct ExactTag : public EffectiveTag { static constexpr ParadigmCode code() { return ParadigmCode::EXACT; } };

struct MetricTag { MetricTag(){}; MetricTag(ValidatedTag){}; static constexpr ParadigmCode code() { return ParadigmCode(7); }  };
struct OrderTag { OrderTag(){}; OrderTag(ValidatedTag){}; static constexpr ParadigmCode code() { return ParadigmCode(3); }  };
struct UpperTag { UpperTag(){}; UpperTag(ValidatedTag){}; static constexpr ParadigmCode code() { return ParadigmCode(2); } };
struct LowerTag { LowerTag(){}; LowerTag(ValidatedTag){}; static constexpr ParadigmCode code() { return ParadigmCode(1); } };
struct ApproximationTag { ApproximationTag(){}; ApproximationTag(ValidatedTag){}; static constexpr ParadigmCode code() { return ParadigmCode(0); } };
using BoundedTag = OrderTag;

class NegatedSierpinskian;
class LowerKleenean;
template<class P> using Equality = P;
template<class P> using LessThan = P;

template<class P1, class P2=P1> struct ParadigmTraits;


template<> struct ParadigmTraits<ExactTag,EffectiveTag> {
    using Weaker = EffectiveTag;
    using Stronger = ExactTag;
};

template<> struct ParadigmTraits<EffectiveTag,EffectiveTag> {
    using Weaker = EffectiveTag;
    using Stronger = EffectiveTag;
    using NextWeaker = ValidatedTag;
};

template<> struct ParadigmTraits<ValidatedTag,ValidatedTag> {
    using Weaker = ValidatedTag;
    using Stronger = ValidatedTag;
    using NextWeaker = ApproximateTag;
};

template<> struct ParadigmTraits<ApproximateTag,ApproximateTag> {
  private:
    typedef ApproximateTag P;
  public:
    using Weaker = ApproximateTag;
    using Stronger = ApproximateTag;
    using NextWeaker = Void;
};

//! \brief A tag meaning that the object provides a positive upper bound for a quantity.
struct ErrorTag { };


//! \brief Inherits from \c TrueType if paradigm \a P1 is weaker than \a P2.
template<class P1, class P2> struct IsWeaker : IsConvertible<P2,P1> { };

//! \brief Inherits from \c TrueType if paradigm \a P1 is stronger than \a P2.
template<class P1, class P2> struct IsStronger : IsWeaker<P2,P1> { };


template<class T> using IsParadigm = IsConvertible<T,ApproximateTag>;

} // namespace Ariadne

#endif /* ARIADNE_PARADIGM_HPP */
