/***************************************************************************
 *            curve.h
 *
 *  Copyright  2007-8  Pieter Collins
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

/*! \file curve.h
 *  \brief A arbitraty curve in Euclidean space.
 */

#ifndef ARIADNE_CURVE_H
#define ARIADNE_CURVE_H

#include <memory>

#include "utility/macros.h"
#include "utility/stlio.h"
#include "function/function.h"
#include "output/graphics_interface.h"

namespace Ariadne {


template<class X> class Vector;

template<class X> class Point;
typedef Point<ExactNumber> ExactPoint;
template<class IVL> class Box;
typedef Box<ExactInterval> ExactBox;
typedef Box<UpperInterval> UpperBox;


// Forward declarations for friends
class CurveInterface;

/*! \brief A curve in Euclidean space
 */
class CurveInterface
{
  public:
    typedef ApproximateFloat ParameterType;
    typedef Point<ApproximateFloat> PointType;
    typedef Vector<ApproximateFloat> TangentVectorType;
  public:
    /*! \brief Destructor. */
    virtual ~CurveInterface() { };
    /*! \brief Return a new dynamically-allocated copy of the curve. */
    virtual CurveInterface* clone() const = 0;
    /*! \brief The dimension of the space the curve lies in. */
    virtual uint dimension() const = 0;
    /*! \brief The smoothness of the curve. */
    virtual ushort smoothness() const = 0;

    /*! \brief The point on the curve at a parameter value. */
    virtual PointType value(const ParameterType& s) const = 0;
    /*! \brief The tangent vector to the curve at a parameter value. */
    virtual TangentVectorType tangent(const ParameterType& s) const = 0;

    /*! \brief Write to an output stream. */
    virtual std::ostream& write(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const CurveInterface& c) {
    return c.write(os); }



/*! \brief A curve in Euclidean space
 */
class Curve
    : public CurveInterface
{
  public:
    typedef CurveInterface::ParameterType ParameterType;
    typedef CurveInterface::PointType PointType;
    typedef CurveInterface::TangentVectorType TangentVectorType;
  public:
    /*! \brief Destructor. */
    virtual ~Curve();
    /*! \brief Constructor. */
    Curve(const EffectiveVectorFunction& f);
    /*! \brief Copy constructor. */
    Curve(const Curve& c);
    /*! \brief Return a new dynamically-allocated copy of the constraint. */
    virtual Curve* clone() const;
    /*! \brief The dimension of the set. */
    virtual uint dimension() const;
    /*! \brief The smoothness of the curve. */
    virtual ushort smoothness() const;

    /*! \brief The value at a point. */
    virtual PointType value(const ParameterType& s) const;
    /*! \brief The tangent at a point. */
    virtual TangentVectorType tangent(const ParameterType& s) const;

    /*! \brief Write to an output stream. */
    virtual std::ostream& write(std::ostream& os) const;
  private:
    EffectiveVectorFunction _function;
};



/*!\brief A line segment in Euclidean space. */
class InterpolatedCurve
    : public DrawableInterface
{
  public:
    typedef CurveInterface::ParameterType ParameterType;
    typedef CurveInterface::PointType PointType;
    typedef CurveInterface::TangentVectorType TangentVectorType;
  public:
    typedef std::map< ParameterType, PointType >::const_iterator const_iterator;
  public:
    /*! \brief Create an empty curve. */
    InterpolatedCurve() : _points() { }
    /*! \brief Create a curve with a single point \a pt at parameter value 0. */
    InterpolatedCurve(const PointType& pt)
        : _points() { this->insert(0,pt); }
    /*! \brief Create a curve with a single point \a pt at parameter value \a s. */
    InterpolatedCurve(ParameterType s, const PointType& pt)
        : _points() { this->insert(s,pt); }
    InterpolatedCurve(const RawFloat& s, const Vector<RawFloat>& pt)
        : _points() { this->insert(s,pt); }
    /*! \brief Create a segment from \a pt0 at parameter value 0 to \a pt1 at parameter value 1. */
    InterpolatedCurve(const PointType& pt0, const PointType& pt1)
        : _points() { this->insert(0,pt0); this->insert(1,pt1); }
    /*! \brief Insert a point with parameter value \a s and spacial value \a pt. */
    void insert(const ParameterType& s, const PointType& pt);
    void insert(const RawFloat& s, const Vector<RawFloat>& pt);

    /*! \brief The number of segments in the curve. */
    size_t size() const { return this->_points.size(); }
    /*! \brief The dimension of the Euclidean space the line segment lies in. */
    uint dimension() const { return this->_points.begin()->second.size(); }
    /*! \brief An iterator to the first point in the curve. */
    const_iterator begin() const { return this->_points.begin(); }
    /*! \brief An iterator to the end point in the curve, NOT the one-past-the-end! */
    const_iterator end() const { return --this->_points.end(); }

    /*! \brief A dynamically-allocated copy. */
    virtual InterpolatedCurve* clone() const;
    /*! \brief Draw on a two-dimensional canvas. */
    virtual void draw(CanvasInterface& c, const Projection2d& p) const;
    /*! \brief A bounding box for the curve. */
    virtual UpperBox bounding_box() const;

    /*! \brief Write to an output stream. */
    virtual std::ostream& write(std::ostream& os) const;

  private:
    friend std::ostream& operator<<(std::ostream&, const InterpolatedCurve&);
  private:
    std::map< ParameterType, PointType > _points;
};

inline
std::ostream& operator<<(std::ostream& os, const InterpolatedCurve& curve) {
    return os << "InterpolatedCurve( size=" << curve.size() << ", points=" << curve._points << " )";
}



} // namespace Ariadne

#endif /* ARIADNE_CURVE_H */