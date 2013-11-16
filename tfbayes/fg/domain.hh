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

#ifndef __TFBAYES_FG_DOMAIN_HH__
#define __TFBAYES_FG_DOMAIN_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cstddef>

#include <boost/icl/interval_set.hpp>

class domain_t {
public:
        domain_t(const boost::icl::interval<double>::type& interval,
                 size_t n) :
                n       (n),
                interval(interval) {
        }
        bool element(const std::vector<double>& x) const {
                assert(x.size() == n);
                for (size_t i = 0; i < x.size(); i++) {
                        if (!boost::icl::contains(interval, x[i])) {
                                return false;
                        }
                }
                return true;
        }
protected:
        size_t n;
        boost::icl::interval<double>::type interval;
};

inline
domain_t real_domain(size_t n) {
        return domain_t(boost::icl::interval<double>::open(
                               -std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity()),
                        n);
}
inline
domain_t positive_domain(size_t n) {
        return domain_t(boost::icl::interval<double>::open(
                                0,
                                std::numeric_limits<double>::infinity()),
                        n);
}

#endif /* __TFBAYES_FG_DOMAIN_HH__ */
