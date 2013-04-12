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
#include <list>
#include <set>
#include <vector>

#include <cassert>
#include <cstdlib>

#include <boost/tuple/tuple.hpp>

#include <phylotree.hh>

class nsplit_t {
public:
        typedef std::vector<size_t> part_t;

        // constructors
        nsplit_t();
        nsplit_t(size_t n, const std::set<size_t>& tmp);

        // methods
        size_t n() const;
        const std::vector<size_t>& part1() const;
        size_t part1(size_t i) const;
        const std::vector<size_t>& part2() const;
        size_t part2(size_t i) const;
        bool null() const;
        // operators
        bool operator<(const nsplit_t& nsplit) const;
        bool operator==(const nsplit_t& nsplit) const;

protected:
        size_t _n;
        part_t _part1;
        part_t _part2;
        bool   _null;
};

class nedge_t : public nsplit_t {
public:
        // constructors
        nedge_t();
        nedge_t(size_t n, std::set<size_t> tmp, double d);
        nedge_t(const nsplit_t& nsplit, double d);

        // methods
        double d() const;
        // operators
        bool operator==(const nedge_t& nedge) const;

protected:
        nsplit_t _nsplit;
        double _d;
};

class common_nedge_t : public nsplit_t {
public:
        // constructors
        common_nedge_t();
        common_nedge_t(const nsplit_t& nsplit, double d1, double d2);

        // methods
        double d1() const;
        double d2() const;

protected:
        nsplit_t _nsplit;
        double _d1;
        double _d2;
};

class nedge_set_t : public std::vector<nedge_t> {
public:
        double length() const;
};

bool compatible(const nsplit_t& s1, const nsplit_t& s2);

class ntree_t {
public:
        // constructors
        ntree_t() { };
        ntree_t(const nedge_set_t& nedge_set,
                const std::vector<double>& leaf_d,
                const std::vector<std::string> leaf_names = std::vector<std::string>());
        ntree_t(const pt_root_t* tree);

        // methods
        pt_root_t* export_tree();
        size_t n() const;
        const nedge_t& find_edge(const nsplit_t& nsplit) const;
        const nedge_set_t& nedge_set() const;
        const nedge_t& nedge_set(size_t i) const;
        const std::vector<double>& leaf_d() const;
        double leaf_d(size_t i) const;
        const std::vector<std::string>& leaf_names() const;
        const std::string& leaf_name(size_t i) const;
        // check splits for compatibility
        bool check_splits() const;

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
        static const std::string _empty_string;
        // empty edge
        static const nedge_t _null_edge;
};

class vertex_cover_t {
public:
        vertex_cover_t(nedge_set_t a, nedge_set_t a_comp,
                       nedge_set_t b, nedge_set_t b_comp,
                       double weight)
                : a(a), a_comp(a_comp), b(b), b_comp(b_comp),
                  weight(weight)
                { }

        // cover for vertices in a
        nedge_set_t a;
        // complement
        nedge_set_t a_comp;
        // cover for vertices in b
        nedge_set_t b;
        // complement
        nedge_set_t b_comp;
        // weight of the vertex cover
        double weight;
};

class incompatibility_graph_t {
public:
         incompatibility_graph_t(const nedge_set_t& a, const nedge_set_t& b);
        ~incompatibility_graph_t();

        vertex_cover_t min_weight_cover(double max_weight = 1.0) const;

        // dimensionality of the constaint matrix
        size_t nrow() const;
        size_t ncol() const;
        // edge sets
        const nedge_set_t& a() const;
        const nedge_set_t& b() const;
        const nedge_t& a(size_t i) const;
        const nedge_t& b(size_t i) const;
        // number of edges in a and b
        size_t na() const;
        size_t nb() const;
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
        nedge_set_t _a;
        nedge_set_t _b;
        size_t _na;
        size_t _nb;

        size_t _nrow;
        size_t _ncol;

        int *_ia, *_ja;
        double *_ar, *_xw;
        bool *_au;
};

class support_pair_t : public std::pair<nedge_set_t, nedge_set_t> {
public:
        support_pair_t()
                : std::pair<nedge_set_t, nedge_set_t>()
                { }
        support_pair_t(const nedge_set_t& s1, const nedge_set_t& s2)
                : std::pair<nedge_set_t, nedge_set_t>(s1, s2)
                { }
        bool operator<(const support_pair_t& s) {
                return first.length()/second.length() < s.first.length()/s.second.length();
        }
};

class npath_t : public std::list<support_pair_t> {
public:
        npath_t() : std::list<support_pair_t>() { }
        npath_t(const support_pair_t& support_pair)
        : std::list<support_pair_t>() {
                insert(begin(), support_pair);
        }
};

class geodesic_t {
public:
        geodesic_t(const ntree_t& t1, const ntree_t& t2);

        ntree_t operator()(const double lambda) const;

        const npath_t& npath() const;
        double length() const;

        const std::list<common_nedge_t>& common_edges() const;
        size_t leaf_n() const;
        const std::vector<std::string>& leaf_names() const;
        const std::vector<double>& t1_leaf_d() const;
        const std::vector<double>& t2_leaf_d() const;
        double t1_leaf_d(size_t i) const;
        double t2_leaf_d(size_t i) const;

protected:
        std::list<common_nedge_t> _common_edges;
        npath_t _npath;
        size_t _leaf_n;
        std::vector<std::string> _leaf_names;
        std::vector<double> _t1_leaf_d;
        std::vector<double> _t2_leaf_d;
};

std::ostream& operator<< (std::ostream& o, const nsplit_t& nsplit);
std::ostream& operator<< (std::ostream& o, const ntree_t& ntree);
std::ostream& operator<< (std::ostream& o, const nedge_t& nedge);
std::ostream& operator<< (std::ostream& o, const nedge_set_t& nedge_set);
std::ostream& operator<< (std::ostream& o, const npath_t& npath);

#endif /* TREESPACE_HH */
