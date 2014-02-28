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

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <glpk.h>

#include <tfbayes/phylotree/treespace.hh>
#include <tfbayes/utility/progress.hh>

using namespace std;

// tools
////////////////////////////////////////////////////////////////////////////////

static
bool empty_intersection(const nsplit_t::part_t& x, const nsplit_t::part_t& y)
{
        nsplit_t::part_t result(x);

        result &= y;

        return !result.any();
}

static
bool is_subset(const nsplit_t::part_t& x, const nsplit_t::part_t& y)
{
        nsplit_t::part_t result(x);

        result &= (~y);

        return !result.any();
}

// nsplit_t
////////////////////////////////////////////////////////////////////////////////

nsplit_t::nsplit_t() : _n(0) { }

nsplit_t::nsplit_t(size_t n, const set<size_t>& tmp)
        : _n(n), _part1(n+1), _part2(n+1)
{
        // switch all bits in part2 on
        _part2.flip();
        // use vectors, which are more comfortable than sets
        vector<size_t> split(tmp.begin(), tmp.end());
        // check arguments
        assert(split.size() > 0);
        assert(split[split.size()-1] <= n);
        // fill part1 and part2
        for (size_t i = 0; i < split.size(); i++) {
                _part1[split[i]] = true;
                _part2[split[i]] = false;
        }
        // leaf zero should always be in part1
        if (_part2[0] == true) {
                _part1.swap(_part2);
        }
}

size_t
nsplit_t::n() const {
        return _n;
}

const nsplit_t::part_t&
nsplit_t::part1() const {
        return _part1;
}

size_t
nsplit_t::part1(size_t i) const {
        return part1()[i];
}

const nsplit_t::part_t&
nsplit_t::part2() const {
        return _part2;
}

size_t
nsplit_t::part2(size_t i) const {
        return part2()[i];
}

bool
nsplit_t::operator==(const nsplit_t& nsplit) const
{
        return part1() == nsplit.part1();
}

bool
nsplit_t::is_ancestor_of(const nsplit_t& nsplit) const
{
        const nsplit_t::part_t& p1 = part2();
        const nsplit_t::part_t& p2 = nsplit.part1();
        const nsplit_t::part_t& p3 = nsplit.part2();

        return !(is_subset(p1, p2) || is_subset(p1, p3));
}

bool compatible(const nsplit_t& s1, const nsplit_t& s2)
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

size_t hash_value(const nsplit_t& nsplit)
{
        return boost::hash_value(nsplit.part1().m_bits);
}

size_t hash_value(const topology_t& topology)
{
        size_t seed = 0;

        for (topology_t::const_iterator it = topology.begin();
             it != topology.end(); it++) {
                // do not use hash combine here, since it depends on
                // the ordering!
                //boost::hash_combine(seed, hash_value(**it));
                seed += hash_value(**it);
        }
        return seed;
}

// nedge_t
////////////////////////////////////////////////////////////////////////////////

nedge_t::nedge_t()
        : nsplit_ptr_t(), _d(0)
{ }

nedge_t::nedge_t(const nedge_t& nedge)
        : nsplit_ptr_t(nedge), _d(nedge.d()), _name(nedge.name())
{ }

nedge_t::nedge_t(size_t n, set<size_t> tmp, double d, string name)
        : nsplit_ptr_t(new nsplit_t(n, tmp)), _d(d), _name(name)
{ }

nedge_t::nedge_t(const nsplit_ptr_t& nsplit_ptr, double d, string name)
        : nsplit_ptr_t(nsplit_ptr), _d(d), _name(name)
{ }

bool
nedge_t::is_ancestor_of(const nedge_t& nedge) const
{
        return operator*().is_ancestor_of(*nedge);
}

const double&
nedge_t::d() const
{
        return _d;
}

double&
nedge_t::d()
{
        return _d;
}

const string&
nedge_t::name() const
{
        return _name;
}

bool
nedge_t::operator==(const nedge_t& nedge) const
{
        return *boost::static_pointer_cast<const nsplit_t>(*this) ==
                *boost::static_pointer_cast<const nsplit_t>(nedge) && d() == nedge.d();
}

// common_nedge_t
////////////////////////////////////////////////////////////////////////////////

common_nedge_t::common_nedge_t()
        : nsplit_ptr_t(), _d1(0), _d2(0)
{ }

common_nedge_t::common_nedge_t(const nsplit_ptr_t& nsplit_ptr, double d1, double d2)
        : nsplit_ptr_t(nsplit_ptr), _d1(d1), _d2(d2)
{ }

double
common_nedge_t::d1() const
{
        return _d1;
}

double
common_nedge_t::d2() const
{
        return _d2;
}

// nedge_set_t
////////////////////////////////////////////////////////////////////////////////

double
nedge_set_t::length() const
{
        double result = 0.0;

        for (nedge_set_t::const_iterator it = begin(); it != end(); it++) {
                result += it->d()*it->d();
        }
        return sqrt(result);
}

pair<nedge_set_t, nedge_set_t>
nedge_set_t::split(nsplit_t split) const
{
        pair<nedge_set_t, nedge_set_t> p;

        for (const_iterator it = begin(); it != end(); it++) {
                const nedge_t& nedge = *it;
                if (split.is_ancestor_of(*nedge)) {
                        p.first.push_back(nedge);
                }
                else {
                        p.second.push_back(nedge);
                }
        }
        return p;
}

// topology_t
////////////////////////////////////////////////////////////////////////////////

topology_t::topology_t(const nedge_set_t& nedge_set)
        : vector<nsplit_ptr_t>() {
        for (nedge_set_t::const_iterator it = nedge_set.begin();
             it != nedge_set.end(); it++) {
                if (it->d() != 0.0) {
                        push_back(*it);
                }
        }
}

bool
topology_t::contains(const nsplit_ptr_t& nsplit_ptr) const
{
        for (const_iterator is = begin(); is != end(); is++) {
                if (*nsplit_ptr == **is) {
                        return true;
                }
        }
        return false;
}

bool
topology_t::operator==(const topology_t& topology) const
{
        if (size() != topology.size()) {
                return false;
        }
        for (const_iterator it = topology.begin(); it != topology.end(); it++) {
                if (!contains(*it)) {
                        return false;
                }
        }
        return true;
}

// nedge_node_t and nedge_root_t
////////////////////////////////////////////////////////////////////////////////

pt_node_t*
convert_leaf_set(
        const vector<double>& leaf_d,
        const vector<string>& leaf_names,
        const nsplit_t::part_t& leaves)
{
        // get the first non-zero bit
        size_t i = leaves.find_first();

        // if there are no leaves, then return NULL
        if (i == leaves.npos) {
                return NULL;
        }

        // start with the first leaf
        pt_leaf_t* current_leaf;
        pt_node_t* current_node = new pt_leaf_t(leaf_d[i], leaf_names[i]);

        // for every other leaf, create a new node and update current_node
        for (i = leaves.find_next(i); i != leaves.npos; i = leaves.find_next(i)) {
                current_leaf = new pt_leaf_t(leaf_d[i], leaf_names[i]);
                current_node = new pt_node_t(0.0, current_node, current_leaf);
        }
        return current_node;
}

nedge_node_t::nedge_node_t()
        : _nedge_set()
{ }

nedge_node_t::nedge_node_t(const map_t& nedge_set)
        : _nedge_set(nedge_set)
{ }

nedge_node_t::nedge_node_t(const nedge_node_t& nedge_node)
        : _nedge_set()
{
        for (map_t::const_iterator it = nedge_node._nedge_set.begin();
             it != nedge_node._nedge_set.end(); it++) {
                _nedge_set[it->first] = it->second->clone();
        }
}

nedge_node_t::~nedge_node_t()
{
        for (map_t::iterator it = _nedge_set.begin();
             it != _nedge_set.end(); it++) {
                delete(it->second);
        }
}

nedge_node_t*
nedge_node_t::clone() const
{
        return new nedge_node_t(*this);
}

nedge_node_t&
nedge_node_t::operator=(const nedge_node_t& nedge_node)
{
        nedge_node_t tmp(nedge_node);
        swap(*this, tmp);
        return *this;
}

bool
nedge_node_t::empty() const
{
        return nedge_set().size() == 0;
}

void
nedge_node_t::propagate(const nedge_t& e)
{
        map_t children;

        // is there an edge which is an ancestor of e?
        for (map_t::iterator it = _nedge_set.begin();
             it != _nedge_set.end(); it++) {
                if (it->first.is_ancestor_of(e)) {
                        it->second->propagate(e);
                        return;
                }
        }
        // is e an ancestor of a subset of edges at this node?
        for (map_t::iterator it = _nedge_set.begin();
             it != _nedge_set.end();) {
                if (e.is_ancestor_of(it->first)) {
                        children[it->first] = it->second;
                        it = _nedge_set.erase(it);
                }
                else {
                        it++;
                }
        }
        _nedge_set[e] = new nedge_node_t(children);
}

nedge_node_t::convert_t
nedge_node_t::convert(
        const vector<double>& leaf_d,
        const vector<string>& leaf_names) const
{
        pt_node_t *tree = NULL, *tmp_node;
        convert_t tmp;
        nsplit_t::part_t tmp_leaves;
        nsplit_t::part_t leaves(leaf_d.size());

        for (map_t::const_iterator it = nedge_set().begin();
             it != nedge_set().end(); it++) {

                if (it->second->empty()) {
                        tmp_node   = convert_leaf_set(leaf_d, leaf_names, it->first->part2());
                        tmp_leaves = it->first->part2();
                        leaves    |= tmp_leaves;
                }
                else {
                        tmp        = it->second->convert(leaf_d, leaf_names);
                        tmp_node   = boost::get<0>(tmp);
                        tmp_leaves = boost::get<1>(tmp);
                        leaves    |= tmp_leaves;
                }
                // create leaves that are missing so far
                if ((it->first->part2() - tmp_leaves).any()) {
                        tmp_node = new pt_node_t(0.0, tmp_node,
                                                 convert_leaf_set(leaf_d, leaf_names, it->first->part2() - tmp_leaves));
                        // add generated leaves
                        leaves  |= it->first->part2() - tmp_leaves;
                }
                // assign edge length
                tmp_node->d = it->first.d();
                // append node to tree
                tree = tree ? new pt_node_t(0.0, tree, tmp_node) : tmp_node;
        }

        return boost::tuples::make_tuple(tree, leaves);
}

nedge_root_t::nedge_root_t(
        const nedge_set_t& nedge_set,
        const vector<double>& leaf_d,
        const vector<string>& leaf_names)
        : nedge_node_t(),
          _leaf_d(leaf_d),
          _leaf_names(leaf_names)
{
        for (size_t i = 0; i < nedge_set.size(); i++) {
                propagate(nedge_set[i]);
        }
}

nedge_root_t::nedge_root_t(const nedge_root_t& nedge_root)
        : nedge_node_t(nedge_root),
          _leaf_d(nedge_root._leaf_d),
          _leaf_names(nedge_root._leaf_names)
{ }

nedge_root_t&
nedge_root_t::operator=(const nedge_root_t& nedge_root)
{
        nedge_root_t tmp(nedge_root);
        swap(*this, tmp);
        return *this;
}

nedge_root_t*
nedge_root_t::clone() const
{
        return new nedge_root_t(*this);
}

pt_root_t
nedge_root_t::convert() const
{
        boost::tuple<pt_node_t*, nsplit_t::part_t> result =
                nedge_node_t::convert(leaf_d(), leaf_names());
        pt_node_t* tree         = boost::get<0>(result);
        nsplit_t::part_t leaves = boost::get<1>(result);

        // leaf 0 should always be missing
        leaves = leaves.flip();
        assert(leaves[0] == 1);
        leaves[0] = 0;

        // create leaves that are missing so far
        if (leaves.any()) {
                if (tree) {
                        // some nodes have been created already
                        tree = new pt_node_t(0.0, tree,
                                             convert_leaf_set(leaf_d(), leaf_names(), leaves));
                }
                else {
                        // star tree topology
                        tree = convert_leaf_set(leaf_d(), leaf_names(), leaves);
                }
        }
        // construct outgroup
        pt_leaf_t* outgroup = new pt_leaf_t(leaf_d(0), leaf_names(0));

        // construct tree
        return pt_root_t(*tree, outgroup);
}

// ntree_t
////////////////////////////////////////////////////////////////////////////////

static set<size_t>
parse_pt_node_t(
        size_t n,
        const pt_node_t& node,
        nedge_set_t& nedge_set,
        vector<double>& leaf_d,
        vector<string>& leaf_names)
{
        set<size_t> s;

        if (node.leaf()) {
                s.insert(node.id);
                leaf_d    [node.id] = node.d;
                leaf_names[node.id] = node.name;
        }
        else {
                set<size_t> s1 = parse_pt_node_t(n, node.left (), nedge_set, leaf_d, leaf_names);
                set<size_t> s2 = parse_pt_node_t(n, node.right(), nedge_set, leaf_d, leaf_names);
                set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(s, s.begin()));
                // insert node only if edge length is greater zero
                if (node.d > 0.0) {
                        nsplit_ptr_t nsplit_ptr(new nsplit_t(n, s));
                        nedge_set.push_back(nedge_t(nsplit_ptr, node.d));
                }
        }
        return s;
}

static void
parse_pt_node_t(
        const pt_root_t& node,
        nedge_set_t& nedge_set,
        vector<double>& leaf_d,
        vector<string>& leaf_names)
{
        assert(node.outgroup());
        // add outgroup
        leaf_d    [node.outgroup()->id] = node.outgroup()->d;
        leaf_names[node.outgroup()->id] = node.outgroup()->name;
        // recursive calls
        parse_pt_node_t(node.n_leaves-1, node.left (), nedge_set, leaf_d, leaf_names);
        parse_pt_node_t(node.n_leaves-1, node.right(), nedge_set, leaf_d, leaf_names);
}

ntree_t::ntree_t(const nedge_set_t& nedge_set,
                 const vector<double>& leaf_d,
                 const vector<string> leaf_names)
        : _n(leaf_d.size()-1),
          _nedge_set(nedge_set),
          _leaf_d(leaf_d),
          _leaf_names(leaf_names)
{
        for (nedge_set_t::const_iterator it = nedge_set.begin(); it != nedge_set.end(); it++) {
                assert(n() == (*it)->n());
        }
        assert(n()+1 == leaf_d.size());
        assert(n()+1 == leaf_names.size());
}

ntree_t::ntree_t(const pt_root_t& tree)
{
        nedge_set_t nedge_set;
        vector<double> leaf_d(tree.n_leaves, 0);
        vector<string> leaf_names(tree.n_leaves, "");

        // parse the tree
        parse_pt_node_t(tree, nedge_set, leaf_d, leaf_names);
        // call constructor
        (*this) = ntree_t(nedge_set, leaf_d, leaf_names);
}

const string ntree_t::_empty_string;
const nedge_t ntree_t::_null_nedge;

const nedge_t&
ntree_t::find_edge(const nsplit_t& nsplit) const
{
        for (nedge_set_t::const_iterator it = nedge_set().begin();
             it != nedge_set().end(); it++) {
                if (*boost::static_pointer_cast<const nsplit_t>(*it) == nsplit) {
                        return *it;
                }
        }
        return _null_nedge;
}

pt_root_t
ntree_t::export_tree() const {
        nedge_root_t nedge_root = nedge_root_t(nedge_set(), leaf_d(), leaf_names());
        return nedge_root.convert();
}

size_t
ntree_t::n() const {
        return _n;
}

const nedge_set_t&
ntree_t::nedge_set() const {
        return _nedge_set;
}

nedge_set_t&
ntree_t::nedge_set() {
        return _nedge_set;
}

const vector<double>&
ntree_t::leaf_d() const {
        return _leaf_d;
}

vector<double>&
ntree_t::leaf_d() {
        return _leaf_d;
}

double
ntree_t::leaf_d(size_t i) const {
        return leaf_d()[i];
}

const vector<string>&
ntree_t::leaf_names() const {
        return _leaf_names;
}

const string&
ntree_t::leaf_name(size_t i) const {
        if (leaf_names().size() > 0) {
                return leaf_names()[i];
        }
        else {
                return _empty_string;
        }
}

bool
ntree_t::compatible(const nsplit_t& nsplit) const
{
        for (nedge_set_t::const_iterator it = nedge_set().begin(); it != nedge_set().end(); it++) {
                if (!::compatible(**it, nsplit)) {
                        return false;
                }
        }
        return true;
}

void
ntree_t::add_edge(const nedge_t& edge)
{
        _nedge_set.push_back(edge);
}

bool
ntree_t::check_splits() const
{
        for (nedge_set_t::const_iterator it = nedge_set().begin(); it != nedge_set().end(); it++) {
                for (nedge_set_t::const_iterator is = nedge_set().begin(); is != it; is++) {
                        if (!::compatible(**it, **is)) {
                                return false;
                        }
                }
        }
        return true;
}

list<common_nedge_t>
ntree_t::common_edges(const ntree_t& tree) const
{
        list<common_nedge_t> result;

        // loop through the set of edges and find common splits
        for (nedge_set_t::const_iterator it = nedge_set().begin();
             it != nedge_set().end(); it++) {
                const nedge_t& edge1 = *it;
                const nedge_t& edge2 = tree.find_edge(*edge1);
                if (edge2) {
                        result.push_back(common_nedge_t(edge1, edge1.d(), edge2.d()));
                }
        }
        return result;
}

void
ntree_t::scale(double lambda)
{
        for (nedge_set_t::iterator it = nedge_set().begin();
             it != nedge_set().end(); it++) {
                it->d() *= lambda;
        }
        for (vector<double>::iterator it = leaf_d().begin();
             it != leaf_d().end(); it++) {
                (*it) *= lambda;
        }
}

topology_t
ntree_t::topology() const
{
        return topology_t(nedge_set());
}

// incompatibility_graph_t
////////////////////////////////////////////////////////////////////////////////

incompatibility_graph_t::incompatibility_graph_t(
        const nedge_set_t& a,
        const nedge_set_t& b)
        : _a(a),
          _b(b),
          _na(a.size()),
          _nb(b.size()),
          _nrow(a.size()*b.size()),
          _ncol(a.size()+b.size())
{
        assert(na() > 0);
        assert(nb() > 0);

        double anorm = pow(a.length(), 2);
        double bnorm = pow(b.length(), 2);

        _ia = (int    *)malloc((1+nrow()*2)*sizeof(int));
        _ja = (int    *)malloc((1+nrow()*2)*sizeof(int));
        _ar = (double *)malloc((1+nrow()*2)*sizeof(double));
        _xw = (double *)malloc((1+ncol()  )*sizeof(double));
        _au = (bool   *)malloc((1+nrow()  )*sizeof(bool));

        // do not use elements indexed by zero
        _ia[0] = 0; _ja[0] = 0;
        _ar[0] = 0; _xw[0] = 0;

        // construct graph
        for (size_t i = 0; i < na(); i++) {
                for (size_t j = 0; j < nb(); j++) {
                        // row index
                        _ia[2*(i*nb() + j) + 1] = i*nb() + j + 1;
                        _ia[2*(i*nb() + j) + 2] = i*nb() + j + 1;
                        // column index
                        _ja[2*(i*nb() + j) + 1] = i + 1;
                        _ja[2*(i*nb() + j) + 2] = j + 1 + na();
                        // value of the constraint matrix
                        if (compatible(*a[i], *b[j])) {
                                _ar[2*(i*nb() + j) + 1] = 0.0;
                                _ar[2*(i*nb() + j) + 2] = 0.0;
                                // edge is not present
                                _au[i*nb() + j] = false;
                        }
                        else {
                                _ar[2*(i*nb() + j) + 1] = 1.0;
                                _ar[2*(i*nb() + j) + 2] = 1.0;
                                // edge is present
                                _au[i*nb() + j] = true;
                        }
                }
        }
        // weights
        for (size_t i = 0; i < na(); i++) {
                // in case the norm is there, there are only edges
                // with length zero, so we have to give them equal
                // weight
                if (anorm == 0) {
                        _xw[i + 1] = 1.0/(double)na();
                }
                else {
                        _xw[i + 1] = pow(a[i].d(), 2)/anorm;
                }
        }
        for (size_t j = 0; j < nb(); j++) {
                if (anorm == 0) {
                        _xw[j + 1 + na()] = 1.0/(double)nb();
                }
                else {
                        _xw[j + 1 + na()] = pow(b[j].d(), 2)/bnorm;
                }
        }
}

incompatibility_graph_t::~incompatibility_graph_t()
{
        free(_ia);
        free(_ja);
        free(_ar);
        free(_xw);
        free(_au);
}

vertex_cover_t
incompatibility_graph_t::min_weight_cover(double max_weight) const
{
        // result
        nedge_set_t ra;
        nedge_set_t ra_comp;
        nedge_set_t rb;
        nedge_set_t rb_comp;
        double weight;
        // integer programming problem
        glp_prob *lp;
        glp_smcp parm;
        lp = glp_create_prob();
        glp_init_smcp(&parm);
        glp_set_obj_dir(lp, GLP_MIN);
        parm.msg_lev = GLP_MSG_OFF;
        // set dimensions
        glp_add_rows(lp, nrow());
        glp_add_cols(lp, ncol());
        // set row and column bounds
        for (size_t i = 0; i < nrow(); i++) {
                if (au(i)) {
                        glp_set_row_bnds(lp, i+1, GLP_LO, 1.0, 0.0);
                }
                else {
                        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
                }
        }
        for (size_t j = 0; j < ncol(); j++) {
                glp_set_col_bnds(lp, j+1, GLP_LO, 0.0, 0.0);
                glp_set_col_kind(lp, j+1, GLP_IV);
        }
        // set weights
        for (size_t k = 0; k < ncol(); k++) {
                glp_set_obj_coef(lp, k+1, xw(k+1));
        }
        // load the constraint matrix
        glp_load_matrix(lp, nrow()*2, ia(), ja(), ar());
        if (glp_simplex(lp, &parm) != 0) {
                cerr << "Warning: The simplex method failed to find a vertex cover!"
                     << endl;
        }
        // save result
        weight = glp_get_obj_val(lp);
        if (weight < max_weight) {
                for (size_t j = 0; j < na(); j++) {
                        glp_get_col_prim(lp, j+1) ? ra.push_back(a(j)) : ra_comp.push_back(a(j));
                }
                for (size_t j = 0; j < nb(); j++) {
                        glp_get_col_prim(lp, na()+j+1) ? rb.push_back(b(j)) : rb_comp.push_back(b(j));
                }
        }
        else {
                ra = a();
                rb = b();
        }
        // free space
        glp_delete_prob(lp);
        // return result
        return vertex_cover_t(ra, ra_comp, rb, rb_comp, weight);
}

const nedge_set_t&
incompatibility_graph_t::a() const
{
        return _a;
}

const nedge_set_t&
incompatibility_graph_t::b() const
{
        return _b;
}

const nedge_t&
incompatibility_graph_t::a(size_t i) const
{
        return _a[i];
}

const nedge_t&
incompatibility_graph_t::b(size_t i) const
{
        return _b[i];
}

size_t
incompatibility_graph_t::na() const
{
        return _na;
}

size_t
incompatibility_graph_t::nb() const
{
        return _nb;
}

size_t
incompatibility_graph_t::nrow() const
{
        return _nrow;
}

size_t
incompatibility_graph_t::ncol() const
{
        return _ncol;
}

int*
incompatibility_graph_t::ia() const
{
        return _ia;
}

int
incompatibility_graph_t::ia(size_t i) const
{
        return _ia[i];
}

int*
incompatibility_graph_t::ja() const
{
        return _ja;
}

int
incompatibility_graph_t::ja(size_t i) const
{
        return _ja[i];
}

double*
incompatibility_graph_t::ar() const
{
        return _ar;
}

double
incompatibility_graph_t::ar(size_t i) const
{
        return _ar[i];
}

double*
incompatibility_graph_t::xw() const
{
        return _xw;
}

double
incompatibility_graph_t::xw(size_t i) const
{
        return _xw[i];
}

bool*
incompatibility_graph_t::au() const
{
        return _au;
}

bool
incompatibility_graph_t::au(size_t i) const
{
        return _au[i];
}

// geodesic_t
////////////////////////////////////////////////////////////////////////////////

geodesic_t::geodesic_t(const ntree_t& __t1, const ntree_t& __t2)
        // start with the cone path
        : //_npath(support_pair_t(t1.nedge_set(), t2.nedge_set())),
        _common_edges(), _npath_list(),
        _leaf_n(__t1.leaf_d().size()),
        _leaf_names(__t1.leaf_names()),
        _t1(__t1),
        _t2(__t2)
{
        // assertions
        assert(t1().leaf_d().size() == t2().leaf_d().size());
        // complement internal edges
        complement_trees();
        // initialize common edges
        _common_edges = t1().common_edges(t2());
        // initialize npath list
        _npath_list = initial_npath_list();
        // apply GTP algorithm to all npaths
        for (list<npath_t>::iterator it = _npath_list.begin(); it != _npath_list.end(); it++) {
                gtp(*it);
        }
}

ntree_t
geodesic_t::operator()(const double lambda) const
{
        // handle sepcial cases
        if (lambda == 0.0) return t1();
        if (lambda == 1.0) return t2();
        // variables
        nedge_set_t nedge_set;
        vector<double> leaf_d(leaf_n(), 0);

        for (list<npath_t>::const_iterator it = npath_list().begin(); it != npath_list().end(); it++) {
                const npath_t& npath(*it);
                npath_t::const_iterator is = npath.begin();
                // fill edge set with edges from B_i
                for (; is != npath.end() &&
                             // slightly different condition than in the paper to prevent
                             // numerical trouble
                             is->first.length()*(1.0-lambda) <= is->second.length()*lambda;
                     is++) {
                        const double length_a = is->first .length();
                        const double length_b = is->second.length();
                        for (nedge_set_t::const_iterator iu = is->second.begin(); iu != is->second.end(); iu++) {
                                const nedge_t& nedge = *iu;
                                const double d = nedge.d()*(lambda*length_b - (1-lambda)*length_a)/length_b;
                                nedge_set.push_back(nedge_t(nedge, d));
                        }
                }
                // fill edge set with edges from A_i
                for (; is != npath.end(); is++) {
                        const double length_a = is->first .length();
                        const double length_b = is->second.length();
                        for (nedge_set_t::const_iterator iu = is->first.begin(); iu != is->first.end(); iu++) {
                                const nedge_t& nedge = *iu;
                                const double d = nedge.d()*((1-lambda)*length_a - lambda*length_b)/length_a;
                                nedge_set.push_back(nedge_t(nedge, d));
                        }
                }
        }
        // add edges that are common to both trees
        for (list<common_nedge_t>::const_iterator it = common_edges().begin();
             it != common_edges().end(); it++) {
                nedge_set.push_back(nedge_t(*it, (1-lambda)*it->d1() + lambda*it->d2()));
        }
        // compute leaf edge lengths
        for (size_t i = 0; i < leaf_d.size(); i++) {
                leaf_d[i] = (1-lambda)*t1_leaf_d(i) + lambda*t2_leaf_d(i);
        }

        return ntree_t(nedge_set, leaf_d, leaf_names());
}

double
geodesic_t::length() const
{
        double result = 0.0;

        // npath distance
        for (list<npath_t>::const_iterator it = npath_list().begin(); it != npath_list().end(); it++) {
                const npath_t& npath(*it);
                for (npath_t::const_iterator is = npath.begin(); is != npath.end(); is++) {
                        result += pow(is->first.length() + is->second.length(), 2);
                }
        }
        // common edge distance
        for (list<common_nedge_t>::const_iterator it = common_edges().begin();
             it != common_edges().end(); it++) {
                result += pow(it->d1() - it->d2(), 2);
        }
        // leaf distance
        for (size_t i = 0; i < leaf_n(); i++) {
                result += pow(t1_leaf_d(i) - t2_leaf_d(i), 2);
        }

        return sqrt(result);
}

const list<common_nedge_t>&
geodesic_t::common_edges() const
{
        return _common_edges;
}

const list<npath_t>&
geodesic_t::npath_list() const
{
        return _npath_list;
}

size_t
geodesic_t::leaf_n() const
{
        return _leaf_n;
}

const vector<string>&
geodesic_t::leaf_names() const
{
        return _leaf_names;
}

const vector<double>&
geodesic_t::t1_leaf_d() const
{
        return t1().leaf_d();
}

const vector<double>&
geodesic_t::t2_leaf_d() const
{
        return t2().leaf_d();
}

double
geodesic_t::t1_leaf_d(size_t i) const
{
        return t1().leaf_d(i);
}

double
geodesic_t::t2_leaf_d(size_t i) const
{
        return t2().leaf_d(i);
}

const ntree_t&
geodesic_t::t1() const
{
        return _t1;
}

const ntree_t&
geodesic_t::t2() const
{
        return _t2;
}

void
geodesic_t::gtp(npath_t& npath)
{
        for (npath_t::iterator it = npath.begin(); it != npath.end();)
        {
                if (it->first.size() == 0 || it->second.size() == 0) {
                        ++it; continue;
                }
                incompatibility_graph_t graph(it->first, it->second);
                vertex_cover_t vc = graph.min_weight_cover();
                if (abs(vc.weight - 1.0) > epsilon && vc.weight < 1.0) {
                        it = npath.erase(it);
                        it = npath.insert(it, support_pair_t(vc.a, vc.b_comp));
                        it = npath.insert(it, support_pair_t(vc.a_comp, vc.b));
                        --it; --it;
                }
                else { ++it; }
        }
        // TODO: sorting is probably not required here!
        npath.sort();
}

const common_nedge_t&
geodesic_t::find_common_edge(const nsplit_t& nsplit) const
{
        for (list<common_nedge_t>::const_iterator it = common_edges().begin();
             it != common_edges().end(); it++) {
                if (**it == nsplit) {
                        return *it;
                }
        }
        return _null_common_nedge;
}

list<npath_t>
geodesic_t::initial_npath_list() const
{
        list<npath_t> result;
        list<nedge_set_t> l1; l1.push_back(nedge_set_t());
        list<nedge_set_t> l2; l2.push_back(nedge_set_t());

        // start with the cone path
        for (nedge_set_t::const_iterator it = t1().nedge_set().begin();
             it != t1().nedge_set().end(); it++) {
                const common_nedge_t& common_nedge = find_common_edge(**it);
                if (!common_nedge) {
                        // this is not a common edge
                        l1.begin()->push_back(*it);
                }
        }
        for (nedge_set_t::const_iterator it = t2().nedge_set().begin();
             it != t2().nedge_set().end(); it++) {
                const common_nedge_t& common_nedge = find_common_edge(**it);
                if (!common_nedge) {
                        // this is not a common edge
                        l2.begin()->push_back(*it);
                }
        }
        // check if there is at least one edge which is not common to
        // both trees
        if (l1.begin()->size() == 0) {
                return result;
        }
        // loop through common edges and split edge sets
        for (list<common_nedge_t>::const_iterator it = common_edges().begin();
             it != common_edges().end(); it++) {
                const common_nedge_t& common_nedge = *it;
                // loop through l1
                for (list<nedge_set_t>::iterator is = l1.begin(); is != l1.end(); is++) {
                        pair<nedge_set_t, nedge_set_t> p = is->split(*common_nedge);
                        if (p.first.size() > 0 && p.second.size() > 0) {
                                is = l1.erase(is);
                                is = l1.insert(is, p.first);
                                is = l1.insert(is, p.second);
                                is--; is--;
                        }
                }
                // loop through l2
                for (list<nedge_set_t>::iterator is = l2.begin(); is != l2.end(); is++) {
                        pair<nedge_set_t, nedge_set_t> p = is->split(*common_nedge);
                        if (p.first.size() > 0 && p.second.size() > 0) {
                                is = l2.erase(is);
                                is = l2.insert(is, p.first);
                                is = l2.insert(is, p.second);
                                is--; is--;
                        }
                }
        }
        for (list<nedge_set_t>::const_iterator it = l1.begin(), is = l2.begin();
             it != l1.end() && is != l2.end(); it++, is++) {
                npath_t npath(support_pair_t(*it, *is));
                result.push_back(npath);
        }
        // loop through l1 and l2 to construct the npath list
        return result;
}

void
geodesic_t::complement_trees()
{
        // one of the trees might have a smaller number of edges, in this
        // case it lies on the boundary between orthants; for the algorithm
        // to work, we need to add those missing edges from the other tree!
        if (t1().nedge_set().size() < t2().nedge_set().size()) {
                for (nedge_set_t::const_iterator it = t2().nedge_set().begin();
                     it != t2().nedge_set().end(); it++) {
                        const nedge_t& nedge = t1().find_edge(**it);
                        if (!nedge && t1().compatible(**it)) {
                                // add edge with length zero
                                _t1.add_edge(nedge_t(*it, 0.0));
                        }
                }
        }
        else if (t2().nedge_set().size() < t1().nedge_set().size()) {
                for (nedge_set_t::const_iterator it = t1().nedge_set().begin();
                     it != t1().nedge_set().end(); it++) {
                        const nedge_t& nedge = t2().find_edge(**it);
                        if (!nedge && t2().compatible(**it)) {
                                // add edge with length zero
                                _t2.add_edge(nedge_t(*it, 0.0));
                        }
                }
        }
}

const common_nedge_t geodesic_t::_null_common_nedge;
const double geodesic_t::epsilon = 0.000001;

// Frechet mean
////////////////////////////////////////////////////////////////////////////////

#include <tfbayes/utility/permutation.hh>

double
frechet_variance(const list<ntree_t>& ntree_list,
                 const vector<double>& weights,
                 const ntree_t& mean)
{
        double result = 0.0;
        double sum    = 0.0;
        // assure that we have as many weights as we have trees
        assert(ntree_list.size() == weights.size());

        list  <ntree_t>::const_iterator it = ntree_list.begin();
        vector<double >::const_iterator is = weights.begin();
        for (; it != ntree_list.end() && is != weights.end(); it++, is++) {
                geodesic_t geodesic(mean, *it);
                result += (*is)*pow(geodesic.length(), 2);
                sum    += *is;
        }
        return result/sum;
}

double
frechet_variance(const list<ntree_t>& ntree_list, const ntree_t& mean)
{
        const vector<double> weights(ntree_list.size(), 1.0);

        return frechet_variance(ntree_list, weights, mean);
}

ntree_t
mean_tree_cyc(const list<ntree_t>& ntree_list, const vector<double>& weights,
              size_t n, const lambda_t& lambda, bool verbose)
{
        // assure that we have as many weights as we have trees
        assert(ntree_list.size() == weights.size());
        // start with the laste element
        ntree_t sk = ntree_list.back();

        for (size_t i = 1, k = 0; k < n; k++) {
                list  <ntree_t>::const_iterator it = ntree_list.begin();
                vector<double >::const_iterator is = weights.begin();
                for (; it != ntree_list.end() && is != weights.end(); it++, is++, i++) {
                        if (verbose && i % 100 == 0) {
                                cerr << progress_t(i/(double)(n*ntree_list.size()));
                        }
                        geodesic_t geodesic(sk, *it);
                        // lambda(k+1) should be constant within one cycle
                        sk = geodesic(2.0*lambda(k+1)*(*is)/(1.0+2.0*lambda(k+1)*(*is)));
                }
        }
        return sk;
}

ntree_t
mean_tree_cyc(const list<ntree_t>& ntree_list, size_t n, const lambda_t& lambda,
              bool verbose)
{
        const vector<double> weights(ntree_list.size(), 1.0);

        return mean_tree_cyc(ntree_list, weights, n, lambda, verbose);
}

ntree_t
median_tree_cyc(const list<ntree_t>& ntree_list, const vector<double>& weights,
                size_t n, const lambda_t& lambda, bool verbose)
{
        // assure that we have as many weights as we have trees
        assert(ntree_list.size() == weights.size());
        // start with the laste element
        ntree_t sk = ntree_list.back();

        for (size_t i = 1, k = 0; k < n; k++) {
                list  <ntree_t>::const_iterator it = ntree_list.begin();
                vector<double >::const_iterator is = weights.begin();
                for (; it != ntree_list.end(); it++, i++) {
                        if (verbose && i % 100 == 0) {
                                cerr << progress_t(i/(double)(n*ntree_list.size()));
                        }
                        geodesic_t geodesic(sk, *it);
                        // lambda(k+1) should be constant within one cycle
                        sk = geodesic(min(1.0, lambda(k+1)*(*is)/geodesic.length()));
                }
        }
        return sk;
}

ntree_t
median_tree_cyc(const list<ntree_t>& ntree_list, size_t n, const lambda_t& lambda,
                bool verbose)
{
        const vector<double> weights(ntree_list.size(), 1.0);

        return median_tree_cyc(ntree_list, weights, n, lambda, verbose);
}

ntree_t
mean_tree_rand(const list<ntree_t>& _ntree_list, const vector<double>& _weights,
               size_t n, boost::random::mt19937& gen, const lambda_t& lambda, bool verbose)
{
        vector<ntree_t> ntree_list(_ntree_list.begin(), _ntree_list.end());
        vector<double > weights(_weights);
        // assure that we have as many weights as we have trees
        assert(ntree_list.size() == weights.size());
        // start with the laste element
        ntree_t sk = ntree_list.back();

        for (size_t i = 1, k = 0; k < n; k++) {
                // shuffle elements in each iteration
                random_permutation_t permutation(ntree_list.size(), gen);
                random_shuffle(ntree_list.begin(), ntree_list.end(), permutation);
                random_shuffle(   weights.begin(),    weights.end(), permutation);

                vector<ntree_t>::const_iterator it = ntree_list.begin();
                vector<double >::const_iterator is = weights.begin();
                for (; it != ntree_list.end() && is != weights.end(); it++, is++, i++) {
                        if (verbose && i % 100 == 0) {
                                cerr << progress_t(i/(double)(n*ntree_list.size()));
                        }
                        geodesic_t geodesic(sk, *it);
                        // in the random version lambda(i+1) changes
                        // in each iteration
                        sk = geodesic(2.0*lambda(i+1)*(*is)/(1.0+2.0*lambda(i+1)*(*is)));
                }
        }
        return sk;
}

ntree_t
mean_tree_rand(const list<ntree_t>& ntree_list, size_t n, boost::random::mt19937& gen, 
               const lambda_t& lambda, bool verbose)
{
        const vector<double> weights(ntree_list.size(), 1.0);

        return mean_tree_rand(ntree_list, weights, n, gen, lambda, verbose);
}

ntree_t
median_tree_rand(const list<ntree_t>& _ntree_list, const vector<double>& _weights,
                 size_t n, boost::random::mt19937& gen, const lambda_t& lambda, bool verbose)
{
        vector<ntree_t> ntree_list(_ntree_list.begin(), _ntree_list.end());
        vector<double > weights(_weights);
        // random generator for shuffling the tree list
        assert(ntree_list.size() == weights.size());
        // start with the laste element
        ntree_t sk = ntree_list.back();

        for (size_t i = 1, k = 0; k < n; k++) {
                // shuffle elements in each iteration
                random_permutation_t permutation(ntree_list.size(), gen);
                random_shuffle(ntree_list.begin(), ntree_list.end(), permutation);
                random_shuffle(   weights.begin(),    weights.end(), permutation);

                vector<ntree_t>::const_iterator it = ntree_list.begin();
                vector<double >::const_iterator is = weights.begin();
                for (; it != ntree_list.end(); it++, i++) {
                        if (verbose && i % 100 == 0) {
                                cerr << progress_t(i/(double)(n*ntree_list.size()));
                        }
                        geodesic_t geodesic(sk, *it);
                        // in the random version lambda(i+1) changes
                        // in each iteration
                        sk = geodesic(min(1.0, lambda(i+1)*(*is)/geodesic.length()));
                }
        }
        return sk;
}

ntree_t
median_tree_rand(const list<ntree_t>& ntree_list, size_t n, boost::random::mt19937& gen,
                 const lambda_t& lambda, bool verbose)
{
        const vector<double> weights(ntree_list.size(), 1.0);

        return median_tree_rand(ntree_list, weights, n, gen, lambda, verbose);
}

ntree_t
mean_same_topology(const std::list<ntree_t>& ntree_list,
                   bool verbose)
{
        size_t n = ntree_list.size();
        // check that the list is not empty
        assert(n > 0);
        // initialize iterator and first mean
        list<ntree_t>::const_iterator it = ntree_list.begin();
        ntree_t mean(*it); it++;

        for (size_t i = 2; it != ntree_list.end(); it++, i++) {
                if (verbose && i % 100 == 0) {
                        cerr << progress_t(i/(double)n);
                }
                mean.scale(i-1);
                mean = geodesic_t(mean, *it)(0.5);
                mean.scale(2.0/(double)i);
        }

        return mean;
}

// consensus trees
////////////////////////////////////////////////////////////////////////////////

ntree_t
majority_consensus(const list<ntree_t>& ntree_list, bool verbose)
{
        assert(ntree_list.size() > 0);

        typedef boost::unordered_map<const nsplit_t, boost::tuple<double, double, nsplit_ptr_t> > nsplit_map_t;

        const size_t n = ntree_list.begin()->n();
        const size_t m = ntree_list.size();
        nsplit_map_t nsplit_map;
        nedge_set_t  nedge_set;
        vector<double> leaf_d(n+1, 0);
        vector<string> leaf_names(ntree_list.begin()->leaf_names());

        // loop through the list of trees an count occurences of splits
        list<ntree_t>::const_iterator it = ntree_list.begin();
        for (size_t i = 1; it != ntree_list.end(); it++, i++) {
                if (verbose && i % 100 == 0) {
                        cerr << progress_t(i/(double)m);
                }
                const ntree_t& ntree(*it);
                for (nedge_set_t::const_iterator is = ntree.nedge_set().begin();
                     is != ntree.nedge_set().end(); is++) {
                        const nedge_t& nedge(*is);
                        if (nsplit_map.find(*nedge) == nsplit_map.end()) {
                                // split not yet present
                                nsplit_map[*nedge] = boost::make_tuple(nedge.d(), 1, nedge);
                        }
                        else {
                                boost::get<0>(nsplit_map[*nedge]) += nedge.d();
                                boost::get<1>(nsplit_map[*nedge]) += 1;
                        }
                }
                // loop over leafs and compute the average branch length
                for (size_t i = 0; i < n+1; i++) {
                        leaf_d[i] += it->leaf_d(i)/(double)m;
                }
        }
        // loop through nsplit_map and find splits that occured in more
        // than half of the trees
        for (nsplit_map_t::const_iterator it = nsplit_map.begin();
             it != nsplit_map.end(); it++) {
                const size_t n = boost::get<1>(it->second);
                const double d = boost::get<0>(it->second)/(double)n;
                if (n > m/2) {
                        nedge_set.push_back(nedge_t(boost::get<2>(it->second), d));
                }
        }
        // construct tree
        return ntree_t(nedge_set, leaf_d, leaf_names);
}

// ostream
////////////////////////////////////////////////////////////////////////////////

ostream&
operator<< (ostream& o, const nsplit_t& nsplit)
{
        o << "{";
        for (nsplit_t::part_t::size_type k = nsplit.part1().find_first();
             k != nsplit.part1().npos;
             k  = nsplit.part1().find_next(k)) {
                if (k != 0) o << ", ";
                o << k;
        }
        o << "} | {";
        for (nsplit_t::part_t::size_type k = nsplit.part2().find_first();
             k != nsplit.part2().npos;
             k = nsplit.part2().find_next(k)) {
                if (k != nsplit.part2().find_first()) o << ", ";
                o << k;
        }
        o << "}";

        return o;
}

ostream&
operator<< (ostream& o, const named_nsplit_t& named_nsplit)
{
        o << "{";
        for (nsplit_t::part_t::size_type k = named_nsplit.part1().find_first();
             k != named_nsplit.part1().npos;
             k  = named_nsplit.part1().find_next(k)) {
                if (k != 0) o << ", ";
                o << named_nsplit.names()[k];
        }
        o << "} | {";
        for (nsplit_t::part_t::size_type k = named_nsplit.part2().find_first();
             k != named_nsplit.part2().npos;
             k  = named_nsplit.part2().find_next(k)) {
                if (k != named_nsplit.part2().find_first()) o << ", ";
                o << named_nsplit.names()[k];
        }
        o << "}";

        return o;
}

ostream&
operator<< (ostream& o, const nedge_t& nedge)
{
        if (nedge.name() != "") {
                o << nedge.name() << ": ";
        }
        o << *nedge << ": "
          << nedge.d();

        return o;
}

ostream&
operator<< (ostream& o, const common_nedge_t& common_nedge)
{
        o << *common_nedge << ": "
          << "(" << common_nedge.d1()
          << "," << common_nedge.d2()
          << ")";

        return o;
}

ostream&
operator<< (ostream& o, const ntree_t& ntree)
{
        o << "internal edges:" << endl
          << ntree.nedge_set()
          << "external edges and leaves:" << endl;
        for (size_t i = 0; i <= ntree.n(); i++) {
                if (ntree.leaf_names().size() > 0) {
                        o << setw(8)
                          << ntree.leaf_name(i)
                          << "(" << i << "): "
                          << ntree.leaf_d(i)
                          << endl;
                }
                else {
                        o << i
                          << ": "
                          << ntree.leaf_d(i)
                          << endl;
                }
        }

        return o;
}

ostream&
operator<< (ostream& o, const nedge_set_t& nedge_set)
{
        for (nedge_set_t::const_iterator it = nedge_set.begin(); it != nedge_set.end(); it++) {
                o << *it << endl;
        }

        return o;
}

ostream& operator<< (ostream& o, const topology_t& topology)
{
        for (topology_t::const_iterator it = topology.begin(); it != topology.end(); it++) {
                o << **it << endl;
        }

        return o;
}

ostream&
operator<< (ostream& o, const npath_t& npath)
{
        for (npath_t::const_iterator it = npath.begin(); it != npath.end(); it++) {
                o << "replacing" << endl << it->first
                  << "with"      << endl << it->second;
        }

        return o;
}

ostream&
operator<< (ostream& o, const list<npath_t>& npath_list)
{
        for (list<npath_t>::const_iterator it = npath_list.begin();
             it != npath_list.end(); it++) {
                o << "npath list element:" << endl
                  << *it;
        }

        return o;
}
