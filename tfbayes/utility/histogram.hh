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

#ifndef __TFBAYES_UTILITY_HISTOGRAM_HH__
#define __TFBAYES_UTILITY_HISTOGRAM_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <ostream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <tfbayes/utility/probability.hh>

template <class input_type = double, class result_type = input_type>
class histogram_t : public std::vector<result_type>
{
        typedef std::vector<result_type> base_t;

        size_t m_n;
        input_type m_min;
        input_type m_max;
        input_type m_width;
        result_type m_total;
        std::vector< input_type> m_x;
        std::vector< input_type> m_counts;
        std::vector<result_type> m_kahan;
public:
        histogram_t()
                : base_t()
                { }
        histogram_t(input_type min, input_type max, size_t n)
                : base_t     (n, 0.0),
                  m_n        (n),
                  m_min      (min),
                  m_max      (max),
                  m_width    ((max-min)/n),
                  m_total    (0.0),
                  m_x        (n, 0.0),
                  m_counts   (n, 0.0),
                  m_kahan    (n, 0.0) {
                for (size_t i = 0; i < n; i++) {
                        m_x[i] = m_width/2.0 + i*m_width;
                }
        }
        template <class T, class S = input_type>
        histogram_t(input_type min, input_type max,
                    const std::vector<T>& y,
                    const std::vector<S>& counts = std::vector<S>())
                : base_t     (y.begin(), y.end()),
                  m_n        (y.size()),
                  m_min      (min),
                  m_max      (max),
                  m_width    ((max-min)/y.size()),
                  m_total    (0.0),
                  m_x        (y.size(), 0.0),
                  m_counts   (y.size(), 0.0),
                  m_kahan    (y.size(), 0.0) {
                for (size_t i = 0; i < y.size(); i++) {
                        // compute midpoints
                        m_x[i] = m_width/2.0 + i*m_width;
                        // sum number of counts
                        m_total += y[i];
                }
                if (std::isnan(m_total)) {
                        throw std::runtime_error("Histogram count overflow!");
                }
                if (counts.size() == y.size()) {
                        std::copy(counts.begin(), counts.end(), m_counts.begin());
                }
        }
        void add(const input_type& x, const result_type& v = 1.0) GCC_ATTRIBUTE_NOAMATH {
                size_t i = which_bin(x);
                assert(i < m_n);
                m_total     += v;
                m_counts[i] += 1.0;
                // kahan summation
                result_type y = v - m_kahan[i];
                result_type t = base_t::operator[](i) + y;
                m_kahan[i] = (t - base_t::operator[](i)) - y;
                base_t::operator[](i) = t;
        }
        size_t which_bin(const input_type& x) const {
                return std::abs(x - m_max) < 1e-8*m_width ? m_n-1 : std::floor((x-m_min)/m_width);
        }
        const size_t& n() const {
                return m_n;
        }
        const input_type& min() const {
                return m_min;
        }
        const input_type& max() const {
                return m_max;
        }
        const input_type& width() const {
                return m_width;
        }
        const result_type& total() const {
                return m_total;
        }
        const std::vector<input_type>& x() const {
                return m_x;
        }
        const std::vector<input_type>& counts() const {
                return m_counts;
        }
        input_type min_counts() const {
                return *std::min_element(m_counts.begin(), m_counts.end());
        }
        friend
        std::ostream& operator<<(std::ostream& o, const histogram_t& hist) {
                input_type min = hist.m_min;
                input_type max = hist.m_min+hist.m_width;
                o << "min max midpoint counts"
                  << std::endl;
                for (size_t i = 0; i < hist.m_n; i++) {
                        o << boost::format("%0.8f %0.8f %0.8f %d")
                                % min % max % hist.m_x[i] % hist[i]
                          << std::endl;
                        min += hist.m_width;
                        max += hist.m_width;
                }
                return o;
        }
private:
        /* serialization */
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & static_cast<base_t&>(*this);
                ar & m_n;
                ar & m_min;
                ar & m_max;
                ar & m_width;
                ar & m_total;
                ar & m_x;
        }
};

template <class input_type, class result_type>
result_type pdf(const histogram_t<input_type, result_type> histogram, const input_type& x) {
        const result_type m = histogram.total();
        const result_type w = histogram.width();
        // no interpolation in these cases
        if (x <= histogram.x().front()) {
                return histogram.front()/(m*w);
        }
        if (x >= histogram.x().back()) {
                return histogram.back()/(m*w);
        }
        // select bins for the interpolation
        size_t j, i = histogram.which_bin(x);
        if (x > histogram.x()[i]) {
                j = i+1;
        }
        else {
                j = i;
                i = i-1;
        }
        // compute interpolation
        const  input_type x1 = histogram.x()[i];
        const result_type y1 = histogram[i];
        const result_type y2 = histogram[j];
        const result_type n  = y1 + (y2 - y1)*static_cast<result_type>((x - x1)/histogram.width());

        return n/(m*w);
}

#endif /* __TFBAYES_UTILITY_HISTOGRAM_HH__ */
