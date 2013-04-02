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

#ifndef TREESPACE_HH
#define TREESPACE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include <cassert>

class nsplit_t {
public:
        nsplit_t(size_t n, std::set<size_t> tmp)
                : _n(n), _part1(tmp.size(), 0), _part2(n-tmp.size(), 0) {
                // use vectors, which are more comfortable
                std::vector<size_t> split(tmp.size(), 0);
                std::copy(tmp.begin(), tmp.end(), split.begin());
                // check arguments
                assert(split.size() > 0);
                assert(split.size() % 2 == 0);
                assert(split[split.size()-1] < n);
                // fill part1 and part2
                for (size_t i = 0, j = 0, k = 0; k < n; k++) {
                        if (split[i] == k) {
                                _part1[i] = k; i++;
                        }
                        else {
                                _part2[j] = k; j++;
                        }
                }
                // leaf zero should always be in part1
                if (_part2[0] == 0) {
                        _part1.swap(_part2);
                }
        }
        size_t n() const {
                return _n;
        }
        const std::vector<size_t>& part1() const {
                return _part1;
        }
        const std::vector<size_t>& part2() const {
                return _part2;
        }
protected:
        size_t _n;
        std::vector<size_t> _part1;
        std::vector<size_t> _part2;
};

bool empty_intersection(const std::vector<size_t>& x, const std::vector<size_t>& y)
{
        // both vectors are assumed to be sorted!
        std::vector<size_t>::const_iterator i = x.begin();
        std::vector<size_t>::const_iterator j = y.begin();
        while (i != x.end() && j != y.end())
        {
                if (*i == *j) {
                        return false;
                }
                else if (*i < *j) {
                        i++;
                }
                else {
                        j++;
                }
        }
        return true;
}

bool compatible(const nsplit_t s1, const nsplit_t s2)
{
        assert(s1.n() == s2.n());
        // s1.part1 and s2.part2 both contain leaf zero, so their
        // intersection is not empty
        if (empty_intersection(s1.part1(), s2.part2())) {
                return true;
        }
        if (empty_intersection(s1.part2(), s2.part1())) {
                return true;
        }
        if (empty_intersection(s1.part2(), s2.part2())) {
                return true;
        }
        return false;
}

std::ostream& operator<< (std::ostream& o, const nsplit_t nsplit);

#endif /* TREESPACE_HH */
