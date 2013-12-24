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

#ifndef _TFBAYES_UTILITY_THREAD_POOL_H_
#define _TFBAYES_UTILITY_THREAD_POOL_H_

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/asio/io_service.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/future.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class thread_pool_t {
public:
        thread_pool_t(size_t n = 1)
                : work(io_service) {
                for (size_t i = 0; i < n; i++) {
                        threads.create_thread(
                                boost::bind(&boost::asio::io_service::run, &io_service)
                                );
                }
        }
        thread_pool_t(const thread_pool_t& thread_pool)
                : work(io_service) {
                for (size_t i = 0; i < thread_pool.threads.size(); i++) {
                        threads.create_thread(
                                boost::bind(&boost::asio::io_service::run, &io_service)
                                );
                }
        }
        virtual ~thread_pool_t() {
                io_service.stop();
                threads.join_all();
        }

        template <typename T>
        boost::unique_future<T> schedule(boost::function<T ()> f) {
                typedef boost::packaged_task<T> task_t;
                boost::shared_ptr<task_t> tmp = boost::make_shared<task_t>(f);
                io_service.post(boost::bind(&task_t::operator(), tmp));
                return tmp->get_future();
        }
protected:
        boost::thread_group threads;
        boost::asio::io_service io_service;
        boost::asio::io_service::work work;
};

#endif /* _TFBAYES_UTILITY_THREAD_POOL_H_ */
