/* Copyright (C) 2015 Philipp Benner
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

#ifndef __TFBAYES_UTILITY_BOOST_RANDOM_SHUFFLE_HH__
#define __TFBAYES_UTILITY_BOOST_RANDOM_SHUFFLE_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>  /* random_shuffle */
#include <functional> /* unary_function */

#include <boost/random/uniform_int_distribution.hpp>

namespace boost { namespace random {

// this class allows to use the boost random number generator
// with the std template library

template<class RandomAccessIterator, class Engine>
void random_shuffle(RandomAccessIterator first, RandomAccessIterator last, Engine& eng)
{
        class std_eng : std::unary_function<unsigned, unsigned> {
        public:
                unsigned operator()(unsigned i) {
                        uniform_int_distribution<> rng(0, i - 1);
                        return rng(m_state);
                }
                std_eng(Engine& state) : m_state(state) {}
        protected:
                Engine &m_state;
        };
        std::random_shuffle(first, last, std_eng(eng));
}


} // namespace random
} // namespace boost


#endif /* __TFBAYES_UTILITY_BOOST_RANDOM_SHUFFLE_HH__ */
