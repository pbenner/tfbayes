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

#ifndef __TFBAYES_UTILITY_RANDOM_HH__
#define __TFBAYES_UTILITY_RANDOM_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <sys/time.h>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/thread.hpp>

class threaded_rng_t : public boost::random::mt19937
{
        typedef boost::random::mt19937 base_t;
public:
        result_type operator()() {
                boost::lock_guard<boost::mutex> guard(mtx);
                return base_t::operator()();
        }
protected:
        boost::mutex mtx;
};

template <typename T>
void seed_rng(T& rng)
{
        // sleep for a millisecond to make sure that we get
        // a unique seed with repeated calls
        boost::this_thread::sleep(boost::posix_time::milliseconds(1));
        // initialize random generator
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

#endif /* __TFBAYES_UTILITY_RANDOM_HH__ */
