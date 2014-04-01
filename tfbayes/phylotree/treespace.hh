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

#ifndef __TFBAYES_PHYLOTREE_TREESPACE_HH__
#define __TFBAYES_PHYLOTREE_TREESPACE_HH__

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

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/utility/clonable.hh>

// nsplit_t forms the basic type for defining edges and ntrees
////////////////////////////////////////////////////////////////////////////////

class nsplit_t {
public:
        typedef boost::dynamic_bitset<> part_t;

        // constructors
        nsplit_t();
        nsplit_t(size_t n, const std::set<size_t>& tmp);

        // methods
        size_t n() const;
        const part_t& part1() const;
        size_t part1(size_t i) const;
        const part_t& part2() const;
        size_t part2(size_t i) const;
        // operators
        bool operator==(const nsplit_t& nsplit) const;


        // ancestor relation where the tree is assumed to be rooted at
        // leaf zero
        bool is_ancestor_of(const nsplit_t& nsplit) const;

protected:
        size_t _n;
        part_t _part1;
        part_t _part2;
};

class named_nsplit_t : public nsplit_t {
public:
        typedef std::vector<std::string> names_t;

        named_nsplit_t(const nsplit_t& nsplit, const names_t& names)
                : nsplit_t(nsplit), _names(names) { }

        const std::vector<std::string> names() const {
                return _names;
        }

        bool operator<(const named_nsplit_t& named_nsplit) const {
                return part1() < named_nsplit.part1();
        }

protected:
        const names_t _names;
};

typedef boost::shared_ptr<const nsplit_t> nsplit_ptr_t;

class nedge_t : public nsplit_ptr_t {
public:
        // constructors
        nedge_t();
        nedge_t(const nedge_t& nedge);
        nedge_t(size_t n, std::set<size_t> tmp, double d, std::string name = "");
        nedge_t(const nsplit_ptr_t& nsplit, double d, std::string name = "");

        // ancestor relation where the tree is assumed to be rooted at
        // leaf zero
        bool is_ancestor_of(const nedge_t& nedge) const;

        // methods
        const double& d() const;
        double& d();
        const std::string& name() const;
        // operators
        bool operator==(const nedge_t& nedge) const;

protected:
        nsplit_t _nsplit;
        double _d;
        std::string _name;
};

// common_nedge_t defines edges that are common to two ntrees
////////////////////////////////////////////////////////////////////////////////

class common_nedge_t : public nsplit_ptr_t {
public:
        // constructors
        common_nedge_t();
        common_nedge_t(const nsplit_ptr_t& nsplit, double d1, double d2);

        // methods
        double d1() const;
        double d2() const;

protected:
        nsplit_t _nsplit;
        double _d1;
        double _d2;
};

// ntrees store their nedges in the nedge_set_t
////////////////////////////////////////////////////////////////////////////////

class nedge_set_t : public std::vector<nedge_t> {
public:
        nedge_set_t()
                : std::vector<nedge_t>() { }
        template <class InputIterator>
        nedge_set_t(InputIterator first, InputIterator last)
                : std::vector<nedge_t>(first, last) { }
        nedge_set_t(const nedge_set_t& nedge_set)
                : std::vector<nedge_t>(nedge_set) { }

        double length() const;
        std::pair<nedge_set_t, nedge_set_t> split(nsplit_t split) const;
};

class topology_t : public std::vector<nsplit_ptr_t> {
public:
        topology_t()
                : std::vector<nsplit_ptr_t>() { }
        topology_t(const nedge_set_t& nedge_set);

        bool contains(const nsplit_ptr_t& nsplit_ptr) const;
        bool operator==(const topology_t& topology) const;
};

bool compatible(const nsplit_t& s1, const nsplit_t& s2);
size_t hash_value(const nsplit_t& nsplit);
size_t hash_value(const topology_t& topology);

// nedge_node_t and nedge_root_t for convertig ntrees to pt_root_t trees
////////////////////////////////////////////////////////////////////////////////

class nedge_node_t : public virtual clonable {
public:
        typedef boost::unordered_map<nedge_t, nedge_node_t*> map_t;
        typedef boost::tuple<pt_node_t*, nsplit_t::part_t> convert_t;

        nedge_node_t();
        nedge_node_t(const map_t& children);
        nedge_node_t(const nedge_node_t& nedge_node);
        virtual ~nedge_node_t();

        virtual nedge_node_t* clone() const;

        virtual nedge_node_t& operator=(const nedge_node_t& nedge_node);

        friend void swap(nedge_node_t& first, nedge_node_t& second) {
                std::swap(first._nedge_set, second._nedge_set);
        }

        virtual bool empty() const;
        virtual void propagate(const nedge_t& e);
        virtual boost::tuple<pt_node_t*, nsplit_t::part_t>
        convert(const std::vector<double>& leaf_d,
                const std::vector<std::string>& leaf_names) const;

        const map_t& nedge_set() const {
                return _nedge_set;
        }

protected:
        // a set of edges where each edge might have a child
        map_t _nedge_set;
};

class nedge_root_t : public nedge_node_t {
public:
        nedge_root_t(const nedge_set_t& nedge_set,
                     const std::vector<double>& leaf_d,
                     const std::vector<std::string>& leaf_names);
        nedge_root_t(const nedge_root_t& nedge_root);
        virtual ~nedge_root_t() { };

        virtual nedge_root_t* clone() const;

        virtual nedge_root_t& operator=(const nedge_root_t& nedge_root);

        friend void swap(nedge_root_t& first, nedge_root_t& second) {
                swap(static_cast<nedge_node_t&>(first),
                     static_cast<nedge_node_t&>(second));
                std::swap(first._leaf_d,     second._leaf_d);
                std::swap(first._leaf_names, second._leaf_names);
        }

        // import convert function so that it is not hidden by our own
        // definition
        using nedge_node_t::convert;
        virtual pt_root_t convert() const;

        const std::vector<double>& leaf_d() const {
                return _leaf_d;
        }
        double leaf_d(size_t i) const {
                return _leaf_d[i];
        }
        const std::vector<std::string>& leaf_names() const {
                return _leaf_names;
        }
        const std::string& leaf_names(size_t i) const {
                return _leaf_names[i];
        }

protected:
        std::vector<double> _leaf_d;
        std::vector<std::string> _leaf_names;
};

// ntree_t definition
////////////////////////////////////////////////////////////////////////////////

class ntree_t {
public:
        typedef std::vector<double> leaf_d_t;
        typedef std::vector<std::string> leaf_names_t;

        // constructors
        ntree_t() : _n(0) { };
        ntree_t(const nedge_set_t& nedge_set,
                const leaf_d_t& leaf_d,
                const leaf_names_t leaf_names = leaf_names_t());
        ntree_t(const pt_root_t& tree);

        // methods
        pt_root_t export_tree() const;
        size_t n() const;
        const nedge_t& find_edge(const nsplit_t& nsplit) const;
        const nedge_set_t& nedge_set() const;
              nedge_set_t& nedge_set();
        const leaf_d_t& leaf_d() const;
              leaf_d_t& leaf_d();
        double leaf_d(size_t i) const;
        const leaf_names_t& leaf_names() const;
        const std::string& leaf_name(size_t i) const;
        bool compatible(const nsplit_t& nsplit) const;
        void add_edge(const nedge_t& edge);
        // check splits for compatibility
        bool check_splits() const;
        // find common edges of two trees, here an edge is common if
        // the split is the same, regardless of the edge length
        std::list<common_nedge_t> common_edges(const ntree_t& tree) const;
        // return the topology of this tree, i.e. the set of splits
        topology_t topology() const;
        // scale the tree by a positive factor
        void scale(double lambda);
        bool null() const { return _n == 0; }

protected:
        size_t _n;
        // n-2 splits
        nedge_set_t _nedge_set;
        // leaf edge lengths, indexed by an integer 0, 1, ..., n
        leaf_d_t _leaf_d;
        // leaf names
        leaf_names_t _leaf_names;
        // empty leaf name
        static const std::string _empty_string;
        // empty edge
        static const nedge_t _null_nedge;
};

// geodesic computations
////////////////////////////////////////////////////////////////////////////////

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
        int* ia() const;
        int* ja() const;
        int  ia(size_t i) const;
        int  ja(size_t i) const;
        // values of the constraint matrix
        double* ar() const;
        double  ar(size_t i) const;
        // vertex weights
        double* xw() const;
        double  xw(size_t i) const;
        // logical values that indicate the presence of an edge
        bool* au() const;
        bool  au(size_t i) const;

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
                const double A1 = first .length();
                const double B1 = second.length();
                const double A2 = s.first .length();
                const double B2 = s.second.length();
                if (B1 == 0.0) {
                        return false;
                }
                if (B2 == 0.0 /* && B1 != 0.0 */) {
                        return true;
                }
                return A1/B1 < A2/B2;
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

        const std::list<npath_t>& npath_list() const;
        double length() const;

        const std::list<common_nedge_t>& common_edges() const;
        size_t leaf_n() const;
        const ntree_t::leaf_names_t& leaf_names() const;
        const ntree_t::leaf_d_t& t1_leaf_d() const;
        const ntree_t::leaf_d_t& t2_leaf_d() const;
        double t1_leaf_d(size_t i) const;
        double t2_leaf_d(size_t i) const;
        const ntree_t& t1() const;
        const ntree_t& t2() const;

protected:
        std::list<npath_t> initial_npath_list() const;
        const common_nedge_t& find_common_edge(const nsplit_t& nsplit) const;
        void gtp(npath_t& npath);
        void complement_trees();

        std::list<common_nedge_t> _common_edges;
        std::list<npath_t> _npath_list;
        size_t _leaf_n;
        ntree_t::leaf_names_t _leaf_names;
        ntree_t _t1;
        ntree_t _t2;
        // empty common edge
        static const common_nedge_t _null_common_nedge;
        // epsilon for comparing double values
        static const double epsilon;
};

// mean and median algorithms
////////////////////////////////////////////////////////////////////////////////

class lambda_t {
public:
        virtual double operator()(size_t k) const = 0;
};

class default_lambda_t : public lambda_t {
public:
        default_lambda_t(double step_size = 1.0)
                : step_size(step_size) { }

        virtual double operator()(size_t k) const {
                return 1.0/(2.0*(1.0 + static_cast<double>(k)/step_size));
        }
        double step_size;
};

double mean_loss(const std::list<ntree_t>& ntree_list, const ntree_t& tree);
double median_loss(const std::list<ntree_t>& ntree_list, const ntree_t& tree);

double frechet_confidence(const std::list<ntree_t>& ntree_list,
                          const ntree_t& mean,
                          double confidence_level);

double frechet_variance(const std::list<ntree_t>& ntree_list,
                        const std::vector<double>& weights,
                        const ntree_t& mean);
double frechet_variance(const std::list<ntree_t>& ntree_list, const ntree_t& mean);

ntree_t mean_tree_cyc(const std::list<ntree_t>& ntree_list, size_t n,
                      const lambda_t& lambda = default_lambda_t(),
                      size_t verbose = false);
ntree_t mean_tree_cyc(const std::list<ntree_t>& ntree_list, const std::vector<double>& weights,
                      size_t n, const lambda_t& lambda = default_lambda_t(),
                      size_t verbose = false);
ntree_t median_tree_cyc(const std::list<ntree_t>& ntree_list, size_t n,
                        const lambda_t& lambda = default_lambda_t(),
                        size_t verbose = false);
ntree_t median_tree_cyc(const std::list<ntree_t>& ntree_list, const std::vector<double>& weights,
                        size_t n, const lambda_t& lambda = default_lambda_t(),
                        size_t verbose = false);

ntree_t mean_tree_rand(const std::list<ntree_t>& ntree_list, size_t n,
                       boost::random::mt19937& gen,
                       const lambda_t& lambda = default_lambda_t(),
                       size_t verbose = false);
ntree_t mean_tree_rand(const std::list<ntree_t>& ntree_list, const std::vector<double>& weights, size_t n,
                       boost::random::mt19937& gen,
                       const lambda_t& lambda = default_lambda_t(),
                       size_t verbose = false);
ntree_t median_tree_rand(const std::list<ntree_t>& ntree_list, size_t n,
                         boost::random::mt19937& gen,
                         const lambda_t& lambda = default_lambda_t(),
                         size_t verbose = false);
ntree_t median_tree_rand(const std::list<ntree_t>& ntree_list, const std::vector<double>& weights, size_t n,
                         boost::random::mt19937& gen,
                         const lambda_t& lambda = default_lambda_t(),
                         size_t verbose = false);

ntree_t mean_same_topology(const std::list<ntree_t>& ntree_list,
                           size_t verbose = false);

// consensus trees
////////////////////////////////////////////////////////////////////////////////

ntree_t majority_consensus(const std::list<ntree_t>& ntree_list, size_t verbose = false);

// i/o operators
////////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< (std::ostream& o, const nsplit_t& nsplit);
std::ostream& operator<< (std::ostream& o, const named_nsplit_t& names_nsplit);
std::ostream& operator<< (std::ostream& o, const ntree_t& ntree);
std::ostream& operator<< (std::ostream& o, const nedge_t& nedge);
std::ostream& operator<< (std::ostream& o, const common_nedge_t& common_nedge);
std::ostream& operator<< (std::ostream& o, const nedge_set_t& nedge_set);
std::ostream& operator<< (std::ostream& o, const topology_t& topology);
std::ostream& operator<< (std::ostream& o, const npath_t& npath);
std::ostream& operator<< (std::ostream& o, const std::list<npath_t>& npath_list);

#endif /* __TFBAYES_PHYLOTREE_TREESPACE_HH__ */
