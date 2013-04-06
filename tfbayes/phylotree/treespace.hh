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
#include <cstdlib>

#include <boost/tuple/tuple.hpp>

#include <phylotree.hh>

class nsplit_t {
public:
        // constructors
        nsplit_t();
        nsplit_t(size_t n, std::set<size_t> tmp);

        // methods
        size_t n() const;
        const std::vector<size_t>& part1() const;
        size_t part1(size_t i) const;
        const std::vector<size_t>& part2() const;
        size_t part2(size_t i) const;

protected:
        size_t _n;
        std::vector<size_t> _part1;
        std::vector<size_t> _part2;
};

class nedge_t : public nsplit_t {
public:
        // constructors
        nedge_t(size_t n, std::set<size_t> tmp, double d);
        nedge_t(const nsplit_t& nsplit, double d);

        double d() const;

protected:
        nsplit_t _nsplit;
        double _d;
};

class nedge_set_t : public std::vector<nedge_t>
{
public:
        double length() const;
};

bool compatible(const nsplit_t& s1, const nsplit_t& s2);

std::ostream& operator<< (std::ostream& o, const nsplit_t nsplit);
std::ostream& operator<< (std::ostream& o, const nedge_t nedge);

class ntree_t {
public:
        // constructors
        ntree_t(const nedge_set_t& nedge_set,
                const std::vector<double>& leaf_d,
                const std::vector<std::string> leaf_names = std::vector<std::string>());

        // methods
        pt_root_t* export_tree();
        size_t n() const;
        const nedge_set_t& nedge_set() const;
        const nedge_t& nedge_set(size_t i) const;
        const std::vector<double>& leaf_d() const;
        double leaf_d(size_t i) const;
        const std::vector<std::string>& leaf_names() const;
        const std::string& leaf_name(size_t i) const;

protected:
        // hidden methods
        boost::tuple<ssize_t, ssize_t, ssize_t> next_splits(const nsplit_t& nsplit, std::vector<bool>& used);
        pt_node_t* export_subtree(std::vector<bool>& used, size_t i);

        size_t _n;
        // n-2 splits
        nedge_set_t _nedge_set;
        // leaf edge lengths,     indexed by an integer 0, 1, ..., n
        std::vector<double> _leaf_d;
        // leaf names
        std::vector<std::string> _leaf_names;
        // empty leaf name
        const std::string _empty_string;
};

class npath_t {
public:
        std::vector<nedge_set_t> A;
        std::vector<nedge_set_t> B;
};

class incompatibility_graph_t {
public:
         incompatibility_graph_t(const nedge_set_t& a, const nedge_set_t& b);
        ~incompatibility_graph_t();

        void min_weight_cover() const;

        // dimensionality of the constaint matrix
        size_t nrow() const;
        size_t ncol() const;
        // row and column indices
        const int* ia() const;
        const int* ja() const;
        const int  ia(size_t i) const;
        const int  ja(size_t i) const;
        // values of the constraint matrix
        const double* ar() const;
        const double  ar(size_t i) const;
        // vertex weights
        const double* xw() const;
        const double  xw(size_t i) const;
        // logical values that indicate the presence of an edge
        const bool* au() const;
        const bool  au(size_t i) const;

protected:
        size_t _nrow;
        size_t _ncol;

        int *_ia, *_ja;
        double *_ar, *_xw;
        bool *_au;
};

#endif /* TREESPACE_HH */
