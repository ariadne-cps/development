/***************************************************************************
 *            concurrency/loggable_smart_thread.hpp
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

/*! \file concurrency/loggable_smart_thread.hpp
 *  \brief A wrapper for smart handling of a thread
 */

#ifndef ARIADNE_LOGGABLE_SMART_THREAD_HPP
#define ARIADNE_LOGGABLE_SMART_THREAD_HPP

#include <utility>
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#include <string>
#include "../output/logging.hpp"

namespace Ariadne {

class LoggableSmartThread {
  public:
    LoggableSmartThread(std::string name, std::function<void(void)> task) {
        _name = name;
        _thread = std::thread([=,this]() {
                    _id = std::this_thread::get_id();
                    _start_promise.set_value();
                    _activate_future.get();
                    if(_active) task();
                });
        _start_future.get();
    }

    std::thread::id id() const {
        return _id;
    }

    std::string name() const {
        return _name;
    }

    void activate()  {
        if (!_active) {
            _active = true;
            Logger::instance().register_thread(*this);
            _activate_promise.set_value();
        }
    }

    ~LoggableSmartThread() {
        if(!_active) _activate_promise.set_value();
        _thread.join();
        Logger::instance().unregister_thread(*this);
    }

  private:
    std::thread::id _id;
    std::string _name;
    std::atomic<bool> _active = false;
    std::thread _thread;
    std::promise<void> _activate_promise;
    std::future<void> _activate_future = _activate_promise.get_future();
    std::promise<void> _start_promise;
    std::future<void> _start_future = _start_promise.get_future();
};

} // namespace Ariadne

#endif // ARIADNE_LOGGABLE_SMART_THREAD_HPP
