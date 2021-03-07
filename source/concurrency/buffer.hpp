/***************************************************************************
 *            concurrency/buffer.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file concurrency/buffer.hpp
 *  \brief A multiple-thread-safe queue usable as a buffer.
 */

#ifndef ARIADNE_BUFFER_HPP
#define ARIADNE_BUFFER_HPP

#include <utility>
#include <mutex>
#include <condition_variable>
#include <queue>
#include "numeric/logical.hpp"

namespace Ariadne {

template<class E> class Buffer
{
  public:
    Buffer(unsigned int cap) : _capacity(cap), _end_consuming(false) { }

    void push(E const& e) {
        std::unique_lock<std::mutex> locker(mu);
        cond.wait(locker, [this](){return _queue.size() < _capacity;});
        _queue.push(e);
        locker.unlock();
        cond.notify_all();
    }
    E pop() {
        std::unique_lock<std::mutex> locker(mu);
        cond.wait(locker, [this](){return _queue.size() > 0 || _end_consuming;});

        if (_end_consuming) throw std::exception();
        E back = _queue.front();
        _queue.pop();
        locker.unlock();
        cond.notify_all();
        return back;
    }

    SizeType size() const {
        return _queue.size();
    }

    void end_consuming() {
        _end_consuming = true;
        cond.notify_all();
    }

    bool ended_consuming() const {
        return _end_consuming;
    }

private:
    std::mutex mu;
    std::condition_variable cond;
    std::queue<E> _queue;
    const unsigned int _capacity;
    std::atomic<bool> _end_consuming;
};

} // namespace Ariadne

#endif // ARIADNE_BUFFER_HPP
