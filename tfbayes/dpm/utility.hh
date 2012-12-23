/* Copyright (C) 2011 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef UTILITY_HH
#define UTILITY_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <vector>
#include <sstream>

static inline
std::vector<std::string> token(const std::string& str, char t) {
        std::string token;
        std::vector<std::string> tokens;
        std::istringstream iss(str);

        while (getline(iss, token, t)) {
                if (token != "") {
                        tokens.push_back(token);
                }
        }
        return tokens;
}

#include <queue>
#include <boost/thread.hpp>

template<typename D>
class save_queue_t
{
public:
        void push(const D& data) {
                boost::mutex::scoped_lock lock(_mutex);
                _queue.push(data);
        }

        bool empty() const {
                boost::mutex::scoped_lock lock(_mutex);
                return _queue.empty();
        }

        D& front() {
                boost::mutex::scoped_lock lock(_mutex);
                return _queue.front();
        }

        D const& front() const {
                boost::mutex::scoped_lock lock(_mutex);
                return _queue.front();
        }

        void pop() {
                boost::mutex::scoped_lock lock(_mutex);
                _queue.pop();
        }
private:
        std::queue<D> _queue;
        mutable boost::mutex _mutex;
};

#endif /* UTILITY_HH */
