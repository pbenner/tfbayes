/* Copyright (C) 2012 Philipp Benner
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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/phylotree/phylotree-gradient-coefficient.hh>

using namespace std;

size_t
hash_value(const pmut_t& pmut)
{
        return (size_t)pmut.node;
}

ostream& operator<< (ostream& o, const mutation_tree_t& tree) {

        if (tree.leaf()) {
                const pmut_t& pmut = tree.pmut();
                if (!!pmut) {
                        if (tree.constant() != 1.0) {
                                o << tree.constant();
                        }
                        if (pmut.mutation) {
                                o << "M("    << pmut.node->name << ")";
                        }
                        else {
                                o << "(1-M(" << pmut.node->name << "))";
                        }
                }
                else {
                        o << tree.constant();
                }
        }
        else {
                if (tree.constant() != 1.0) {
                        o << tree.constant();
                }
                o << "(" << tree.left() << ")";
                if (tree.mult()) {
                        o << "*";
                }
                if (tree.sum()) {
                        o << "+";
                }
                o << "(" << tree.right() << ")";
        }

        return o;
}
