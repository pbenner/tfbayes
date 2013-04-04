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

using namespace std;

/* tools
 *****************************************************************************/

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

vector<size_t> difference(const vector<size_t>& x, const vector<size_t>& y)
{
        vector<size_t> result;
        // both vectors are assumed to be sorted!
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

/* ntree_t
 *****************************************************************************/

ntree_t::ntree_t(const vector<nsplit_t>& splits,
                 const vector<double>& int_d,
                 const vector<double>& leaf_d,
                 const vector<string> leaf_names)
        : _n(splits.size()+2),
          _splits(splits),
          _int_d(int_d),
          _leaf_d(leaf_d),
          _leaf_names(leaf_names) {
        // check that there is at least one split
        assert(splits.size() > 0);
        // we need n-2 splits to fully specify a tree
        assert(splits.size() == splits.begin()->n()-2);
        assert(leaf_names.size() == 0 || leaf_names.size() == n()+1);
        // check splits for compatibility
        for (vector<nsplit_t>::const_iterator it = splits.begin(); it != splits.end(); it++) {
                for (vector<nsplit_t>::const_iterator is = splits.begin(); is != it; is++) {
                        if (!compatible(*it, *is)) {
                                cerr << "Invalid set of splits."
                                          << endl;
                                exit(EXIT_FAILURE);
                        }
                }
        }
}

boost::tuple<ssize_t, ssize_t, ssize_t>
ntree_t::next_splits(const nsplit_t& split, vector<bool>& used) {
        // return value:
        // (int edge, int edge, -1) OR (int edge, -1, leaf edge)
        //
        // check first if we need only one internal edge
        for (size_t i = 0; i < n()-2; i++) {
                if (used[i]) continue;
                if (splits(i).part1().size() == split.part1().size() + 1) {
                        vector<size_t> tmp = difference(splits(i).part1(), split.part1());
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
                        vector<size_t> tmp = intersection(splits(i).part1(), splits(j).part1());
                        if (equal(tmp.begin(), tmp.end(), split.part1().begin())) {
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
        if (splits(i).part2().size() == 2) {
                left_tree  = new pt_leaf_t(-1, leaf_d(splits(i).part2(0)), leaf_name(splits(i).part2(0)));
                right_tree = new pt_leaf_t(-1, leaf_d(splits(i).part2(1)), leaf_name(splits(i).part2(1)));
        }
        // otherwise we need to do a recursive call
        else {
                // find the next internal edge(s)
                boost::tuple<ssize_t, ssize_t, ssize_t> ns = next_splits(splits(i), used);
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
        return new pt_node_t(-1, int_d(i), left_tree, right_tree);
}

size_t
ntree_t::n() const {
        return _n;
}

const vector<nsplit_t>&
ntree_t::splits() const {
        return _splits;
}

const nsplit_t&
ntree_t::splits(size_t i) const {
        return splits()[i];
}

const vector<double>&
ntree_t::int_d() const {
        return _int_d;
}

double
ntree_t::int_d(size_t i) const {
        return int_d()[i];
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
