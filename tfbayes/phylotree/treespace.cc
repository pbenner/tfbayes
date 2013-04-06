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

#include <treespace.hh>

#include <cmath>
#include <cstdlib>

#include <glpk.h>

using namespace std;

/* tools
 *****************************************************************************/

static
bool empty_intersection(const vector<size_t>& x, const vector<size_t>& y)
{
        // both vectors are assumed to be sorted!
        vector<size_t>::const_iterator i = x.begin();
        vector<size_t>::const_iterator j = y.begin();
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

static
vector<size_t> intersection(const vector<size_t>& x, const vector<size_t>& y)
{
        vector<size_t> result;
        // both vectors are assumed to be sorted!
        vector<size_t>::const_iterator i = x.begin();
        vector<size_t>::const_iterator j = y.begin();
        while (i != x.end() && j != y.end())
        {
                if (*i == *j) {
                        result.push_back(*i);
                        i++; j++;
                }
                else if (*i < *j) {
                        i++;
                }
                else {
                        j++;
                }
        }
        return result;
}

static
vector<size_t> difference(const vector<size_t>& x, const vector<size_t>& y)
{
        vector<size_t> result;

        set_difference(x.begin(), x.end(), y.begin(), y.end(),
                       back_inserter(result));

        return result;
}

/* nsplit_t
 *****************************************************************************/

nsplit_t::nsplit_t() : _n(0) { }

nsplit_t::nsplit_t(size_t n, set<size_t> tmp)
        : _n(n), _part1(tmp.size(), 0), _part2(n-tmp.size()+1, 0) {
        // use vectors, which are more comfortable
        vector<size_t> split(tmp.size(), 0);
        copy(tmp.begin(), tmp.end(), split.begin());
        // check arguments
        assert(split.size() > 0);
        assert(split[split.size()-1] <= n);
        // fill part1 and part2
        for (size_t i = 0, j = 0, k = 0; k <= n; k++) {
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

size_t
nsplit_t::n() const {
        return _n;
}

const vector<size_t>&
nsplit_t::part1() const {
        return _part1;
}

size_t
nsplit_t::part1(size_t i) const {
        return part1()[i];
}

const vector<size_t>&
nsplit_t::part2() const {
        return _part2;
}

size_t
nsplit_t::part2(size_t i) const {
        return part2()[i];
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

ostream&
operator<< (ostream& o, const nsplit_t nsplit)
{
        o << "{";
        for (size_t i = 0; i < nsplit.part1().size(); i++) {
                if (i < nsplit.part1().size()-1) {
                        o << nsplit.part1(i)
                          << ", ";
                }
                else {
                        o << nsplit.part1(i)
                          << "} | {";
                }
        }
        for (size_t i = 0; i < nsplit.part2().size(); i++) {
                if (i < nsplit.part2().size()-1) {
                        o << nsplit.part2()[i]
                          << ", ";
                }
                else {
                        o << nsplit.part2()[i]
                          << "}";
                }
        }
        return o;
}

ostream&
operator<< (ostream& o, const nedge_t nedge)
{
        o << nedge.d() << ": "
          << static_cast<nsplit_t>(nedge);

        return o;
}

/* nedge_t
 *****************************************************************************/

nedge_t::nedge_t(size_t n, std::set<size_t> tmp, double d)
        : nsplit_t(n, tmp), _d(d)
{ }

nedge_t::nedge_t(const nsplit_t& nsplit, double d)
        : nsplit_t(nsplit), _d(d)
{ }

double
nedge_t::d() const
{
        return _d;
}

/* nedge_set_t
 *****************************************************************************/

double
nedge_set_t::length() const
{
        double result = 0.0;

        for (nedge_set_t::const_iterator it = begin(); it != end(); it++) {
                result += it->d()*it->d();
        }
        return sqrt(result);
}

/* ntree_t
 *****************************************************************************/

ntree_t::ntree_t(const nedge_set_t& nedge_set,
                 const vector<double>& leaf_d,
                 const vector<string> leaf_names)
        : _n(nedge_set.size()+2),
          _nedge_set(nedge_set),
          _leaf_d(leaf_d),
          _leaf_names(leaf_names) {
        // check that there is at least one split
        assert(nedge_set.size() > 0);
        // we need n-2 internal edges to fully specify a tree
        assert(nedge_set.size() == nedge_set.begin()->n()-2);
        assert(leaf_names.size() == 0 || leaf_names.size() == n()+1);
        // check splits for compatibility
        for (nedge_set_t::const_iterator it = nedge_set.begin(); it != nedge_set.end(); it++) {
                for (nedge_set_t::const_iterator is = nedge_set.begin(); is != it; is++) {
                        if (!compatible(*it, *is)) {
                                cerr << "Invalid set of internal edges."
                                          << endl;
                                exit(EXIT_FAILURE);
                        }
                }
        }
}

boost::tuple<ssize_t, ssize_t, ssize_t>
ntree_t::next_splits(const nsplit_t& nsplit, vector<bool>& used) {
        // return value:
        // (int edge, int edge, -1) OR (int edge, -1, leaf edge)
        //
        // check first if we need only one internal edge
        for (size_t i = 0; i < n()-2; i++) {
                if (used[i]) continue;
                if (nedge_set(i).part1().size() == nsplit.part1().size() + 1) {
                        vector<size_t> tmp = difference(nedge_set(i).part1(), nsplit.part1());
                        if (tmp.size() == 1) {
                                // mark edge used
                                used[i] = true;
                                // return index of that edge
                                return boost::make_tuple(i, -1, tmp[0]);
                        }
                }
        }
        // now search for the two next internal edges
        for (size_t i = 0; i < n()-2; i++) {
                if (used[i]) continue;
                for (size_t j = 0; j < i; j++) {
                        if (used[j]) continue;
                        vector<size_t> tmp = intersection(nedge_set(i).part1(), nedge_set(j).part1());
                        if (equal(tmp.begin(), tmp.end(), nsplit.part1().begin())) {
                                // mark both edges used
                                used[i] = true;
                                used[j] = true;
                                // return indices of both edges
                                return boost::make_tuple(i, j, -1);
                        }
                }
        }
        // we should never arrive here
        return boost::make_tuple(-1, -1, -1);
}

pt_root_t*
ntree_t::export_tree() {
        pt_node_t* left_tree;
        pt_node_t* right_tree;
        // the vecor used stores which edges were already used
        // (this should optimize performance)
        vector<bool>   used(n()-2, false);
        // list of leafs for a given subtree
        vector<size_t> leafs(n(), 0);
        for (size_t i = 0; i <= n(); i++) {
                leafs[i] = i;
        }
        // start with leaf zero to build the tree
        set<size_t> tmp; tmp.insert(0);
        boost::tuple<ssize_t, ssize_t, ssize_t> ns = next_splits(nsplit_t(n(), tmp), used);
        if (boost::get<1>(ns) != -1) {
                left_tree  = export_subtree(used, boost::get<0>(ns));
                right_tree = export_subtree(used, boost::get<1>(ns));
        }
        else {
                left_tree  = export_subtree(used, boost::get<0>(ns));
                right_tree = new pt_leaf_t(-1, leaf_d(boost::get<2>(ns)), leaf_name(boost::get<2>(ns)));
        }
        // return resulting tree
        return new pt_root_t(-1, left_tree, right_tree);
}

pt_node_t*
ntree_t::export_subtree(vector<bool>& used, size_t i) {
        pt_node_t* left_tree;
        pt_node_t* right_tree;
        // check if there are only leafs following
        if (nedge_set(i).part2().size() == 2) {
                left_tree  = new pt_leaf_t(-1, leaf_d(nedge_set(i).part2(0)), leaf_name(nedge_set(i).part2(0)));
                right_tree = new pt_leaf_t(-1, leaf_d(nedge_set(i).part2(1)), leaf_name(nedge_set(i).part2(1)));
        }
        // otherwise we need to do a recursive call
        else {
                // find the next internal edge(s)
                boost::tuple<ssize_t, ssize_t, ssize_t> ns = next_splits(nedge_set(i), used);
                // check if there is one or two edges folliwing
                if (boost::get<1>(ns) != -1) {
                        left_tree  = export_subtree(used, boost::get<0>(ns));
                        right_tree = export_subtree(used, boost::get<1>(ns));
                }
                else {
                        left_tree  = export_subtree(used, boost::get<0>(ns));
                        right_tree = new pt_leaf_t(-1, leaf_d(boost::get<2>(ns)), leaf_name(boost::get<2>(ns)));
                }
        }
        return new pt_node_t(-1, nedge_set(i).d(), left_tree, right_tree);
}

size_t
ntree_t::n() const {
        return _n;
}

const nedge_set_t&
ntree_t::nedge_set() const {
        return _nedge_set;
}

const nedge_t&
ntree_t::nedge_set(size_t i) const {
        return nedge_set()[i];
}

const vector<double>&
ntree_t::leaf_d() const {
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

/* incompatibility_graph_t
 *****************************************************************************/

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

        _ia = (int    *)malloc(1+nrow()*2*sizeof(int));
        _ja = (int    *)malloc(1+nrow()*2*sizeof(int));
        _ar = (double *)malloc(1+nrow()*2*sizeof(double));
        _xw = (double *)malloc(1+ncol()*  sizeof(double));
        _au = (bool   *)malloc(1+nrow()*  sizeof(bool));

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
                        if (compatible(a[i], b[j])) {
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
                _xw[i + 1] = pow(a[i].d(), 2)/anorm;
        }
        for (size_t j = 0; j < nb(); j++) {
                _xw[j + 1 + na()] = pow(b[j].d(), 2)/bnorm;
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
        glp_simplex(lp, &parm);
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

const int*
incompatibility_graph_t::ia() const
{
        return _ia;
}

const int
incompatibility_graph_t::ia(size_t i) const
{
        return _ia[i];
}

const int*
incompatibility_graph_t::ja() const
{
        return _ja;
}

const int
incompatibility_graph_t::ja(size_t i) const
{
        return _ja[i];
}

const double*
incompatibility_graph_t::ar() const
{
        return _ar;
}

const double
incompatibility_graph_t::ar(size_t i) const
{
        return _ar[i];
}

const double*
incompatibility_graph_t::xw() const
{
        return _xw;
}

const double
incompatibility_graph_t::xw(size_t i) const
{
        return _xw[i];
}

const bool*
incompatibility_graph_t::au() const
{
        return _au;
}

const bool
incompatibility_graph_t::au(size_t i) const
{
        return _au[i];
}