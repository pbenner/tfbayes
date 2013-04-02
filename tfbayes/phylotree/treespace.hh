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
#include <vector>

#include <cassert>

class nsplit_t {
public:
        nsplit_t(size_t n, std::vector<size_t> split)
                : _n(n), _part1(split.size(), 0), _part2(n-split.size(), 0) {
                // sort split so we can easily fill part1 and part2
                std::sort(split.begin(), split.end());
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

std::ostream& operator<< (std::ostream& o, const nsplit_t nsplit);

#endif /* TREESPACE_HH */
