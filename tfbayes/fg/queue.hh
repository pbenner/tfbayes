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

#include <set>

#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>

#include <tfbayes/fg/node-types.hh>

class partial_result_t : std::map<const fg_node_i*, double> {
public:
        typedef std::map<const fg_node_i*, double> base_t;

        partial_result_t() :
                sum(0.0) {
        }
        void update(const fg_node_i* node, double new_value) {
                mtx.lock();
                double old_value = base_t::operator[](node);
                sum -= old_value;
                sum += new_value;
                base_t::operator[](node) = new_value;
                mtx.unlock();
        }
        double operator()() const {
                return sum;
        }
protected:
        double sum;
        boost::mutex mtx;
};

template <typename T>
class node_queue_t : std::set<T*> {
public:
        typedef std::set<T*> base_t;

        using base_t::empty;
        void push(T* node) {
                mtx.lock();
                base_t::insert(node);
                mtx.unlock();
        }
        T* pop() {
                T* node = NULL;
                mtx.lock();
                typename base_t::iterator it = base_t::begin();
                if (it != base_t::end()) {
                        node = *it;
                        base_t::erase(node);
                }
                mtx.unlock();
                return node;
        }
private:
        boost::mutex mtx;
};

class fg_queue_t {
public:
        typedef std::set<fg_node_i*> base_t;

        fg_queue_t(size_t threads) :
                barrier(threads),
                // variable nodes should go first
                which  (false),
                first1 (false),
                first2 (false),
                done   (false) {
        }
        void set_limit(boost::optional<size_t> n = boost::optional<size_t>()) {
                limit = n;
        }
        void push_factor(factor_node_i* node) {
                mtx.lock();
                factor_queue.push(node);
                mtx.unlock();
        }
        void push_variable(variable_node_i* node) {
                mtx.lock();
                variable_queue.push(node);
                mtx.unlock();
        }
        fg_node_i* pop() {
                fg_node_i* node = NULL;
                while (node == NULL) {
                        mtx.lock();
                        if (which && factor_queue.empty()) {
                                first1 = true;
                                mtx.unlock();
                                // barrier (wait for all other threads to
                                // arrive here)
                                barrier.wait();
                                mtx.lock();
                                // is this thread the first to exit
                                // the barrier?
                                if (first1) {
                                        first1 = false;
                                        // check if all jobs are done
                                        if (variable_queue.empty()) {
                                                done = true;
                                        }
                                        if (limit && *limit == 0) {
                                                done = true;
                                        }
                                        debug("################################################################################"
                                              << std::endl);
                                        // switch to other queue
                                        which = false;
                                        // save free energy
                                        debug(boost::format("summed free energy is: %d\n")
                                              % partial_result());
                                        history.push_back(partial_result());
                                }
                        }
                        if (!which && variable_queue.empty()) {
                                first2 = true;
                                mtx.unlock();
                                // barrier (wait for all other threads to
                                // arrive here)
                                barrier.wait();
                                mtx.lock();
                                // is this thread the first to exit
                                // the barrier?
                                if (first2) {
                                        first2 = false;
                                        // check if all jobs are done
                                        if (factor_queue.empty()) {
                                                done = true;
                                        }
                                        if (limit && *limit == 0) {
                                                done = true;
                                        }
                                        debug("################################################################################"
                                              << std::endl);
                                        // switch to other queue
                                        which = true;
                                }
                        }
                        if (done) {
                                mtx.unlock();
                                return NULL;
                        }
                        // get job
                        node = which
                                ? static_cast<fg_node_i*>(factor_queue.pop())
                                : static_cast<fg_node_i*>(variable_queue.pop());
                        if (node != NULL) {
                                if (limit && *limit > 0) {
                                        *limit -= 1;
                                }
                        }
                        mtx.unlock();
                }
                return node;
        }
        void save_result(fg_node_i* job, double result) {
                partial_result.update(job, result);
        }
        // record free energy
        std::vector<double> history;
protected:
        // thread barrier
        boost::barrier barrier;
        // from which queue to receive jobs
        bool which;
        // indicators for whether a thread leaves the barrier first
        bool first1;
        bool first2;
        // indicates whether all jobs have been processed
        bool done;
        // maximum number of jobs to process
        boost::optional<size_t> limit;
        // partial results
        partial_result_t partial_result;
        // the two queues
        node_queue_t<  factor_node_i> factor_queue;
        node_queue_t<variable_node_i> variable_queue;
        boost::mutex mtx;
};

#endif /* __TFBAYES_FG_QUEUE_HH__ */
