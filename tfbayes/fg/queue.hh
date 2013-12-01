/* Copyright (C) 2013 Philipp Benner
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

#ifndef __TFBAYES_FG_QUEUE_HH__
#define __TFBAYES_FG_QUEUE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <map>
#include <queue>
#include <set>

#include <tfbayes/fg/node-types.hh>

template <typename T>
class cache_t : std::map<const T*, double> {
public:
        typedef std::map<const T*, double> base_t;

        cache_t() :
                sum(0.0) {
        }
        void update(const T* node, double new_value) {
                double old_value = base_t::operator[](node);
                sum -= old_value;
                sum += new_value;
                base_t::operator[](node) = new_value;
        }
        double operator()() const {
                // save free energy
                debug(boost::format("summed free energy is: %d\n") % sum);
                return sum;
        }
protected:
        double sum;
};

template <typename T>
class fg_queue_t : std::queue<T*> {
public:
        typedef std::queue<T*> base_t;

        void set_limit(boost::optional<size_t> n = boost::optional<size_t>()) {
                limit = n;
        }
        void push(T* node) {
                if (set.find(node) == set.end()) {
                        base_t::push(node);
                        set.insert(node);
                }
        }
        T* pop() {
                T* node = NULL;
                if ((!limit || *limit != 0) && !base_t::empty()) {
                        node = base_t::front();
                               base_t::pop();
                        set.erase(node);
                }
                // update limit
                if (limit && *limit > 0) {
                        *limit -= 1;
                }
                return node;
        }
protected:
        // assure that every job is queued only once
        std::set<T*> set;
        // maximum number of jobs to process
        boost::optional<size_t> limit;
};

#endif /* __TFBAYES_FG_QUEUE_HH__ */
