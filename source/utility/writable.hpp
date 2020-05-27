/***************************************************************************
 *            utility/writable.hpp
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

/*! \file utility/writable.hpp
 *  \brief
 */



#ifndef ARIADNE_WRITABLE_HPP
#define ARIADNE_WRITABLE_HPP

#include "typedefs.hpp"
#include "metaprogramming.hpp"

namespace Ariadne {

/************ WritableInterface **********************************************/

class WritableInterface {
  public:
    virtual ~WritableInterface() = default;
    friend OutputStream& operator<<(OutputStream& os, const WritableInterface& w);
  public:
    inline OutputStream& write(OutputStream& os) const { return this->_write(os); }
  protected:
  public:
    virtual OutputStream& _write(OutputStream&) const = 0;
};
inline OutputStream& operator<<(OutputStream& os, const WritableInterface& w) { w._write(os); return os; }

template<class T, class = decltype(declval<T>()._write(declval<OutputStream>()))> True has_write(int);
template<class T> False has_write(...);
template<class T, class = Fallback> struct IsWritable : decltype(has_write<T>(1)) { };

template<class T> requires IsWritable<T>::value OutputStream& operator<<(OutputStream& os, const T& t) {
    return t._write(os);
}

template<class T> class Writer;
template<class T> class WriterInterface;
template<class T, class W=WriterInterface<T>> class WritableTemporary;

template<class T> class WriterInterface {
  public:
    virtual ~WriterInterface() = default;
    inline WritableTemporary<T> operator() (T const& t) const;
  private:
    virtual OutputStream& _write(OutputStream& os, T const& t) const = 0;
    friend OutputStream& operator<<(OutputStream& os, WritableTemporary<T> const& t);
};

template<class T> class Handle;
template<class T> class Writer : public Handle<WriterInterface<T>> {
  public:
    using Handle<WriterInterface<T>>::Handle;
    Writer<T>(Handle<WriterInterface<T>> wh) : Handle<WriterInterface<T>>(wh) { }
    inline WritableTemporary<T> operator() (T const& t) const;
};


template<class T, class W> class WritableTemporary {
  private:
    W const& _w; T const& _t;
    WritableTemporary(W const& w, T const& t) : _w(w), _t(t) { }
    friend W; friend class Writer<T>;
  public:
    friend OutputStream& operator<<(OutputStream& os, WritableTemporary<T,W> const& wt) { return wt._w._write(os,wt._t); }
};

template<class T> WritableTemporary<T> WriterInterface<T>::operator() (T const& t) const {
    return WritableTemporary<T>(*this,t); }
template<class T> WritableTemporary<T> Writer<T>::operator() (T const& t) const {
    return WritableTemporary<T>(*this->_ptr,t); }


template<class T> class RepresentationWriter;

} // namespace Ariadne

#endif
