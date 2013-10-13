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

#ifndef PHYLOTREE_GENERATE_OBSERVATIONS_HH
#define PHYLOTREE_GENERATE_OBSERVATIONS_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cstdlib>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/utility/statistics.hh>

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
void pt_generate_observations(const pt_node_t& node,
                              const std::vector<double>& stationary,
                              const CODE_TYPE parent_nucleotide,
                              std::vector<CODE_TYPE>& observations)
{
        double r = (double)rand()/RAND_MAX;
        CODE_TYPE nucleotide;
        if (r <= node.mutation_probability()) {
                nucleotide = categorial_sample<CODE_TYPE, ALPHABET_SIZE>(stationary);
        }
        else {
                nucleotide = parent_nucleotide;
        }
        if (node.leaf()) {
                /* save nucleotide in our observations */
                observations[node.id] = nucleotide;
        }
        /* traverse the treee */
        else {
                pt_generate_observations<CODE_TYPE, ALPHABET_SIZE>(
                        node.left (), stationary, nucleotide, observations);
                pt_generate_observations<CODE_TYPE, ALPHABET_SIZE>(
                        node.right(), stationary, nucleotide, observations);
        }
}

template <typename CODE_TYPE, size_t ALPHABET_SIZE>
std::vector<CODE_TYPE>
pt_generate_observations(const pt_root_t& tree, const std::vector<double>& stationary)
{
        std::vector<CODE_TYPE> observations(tree.n_leaves, -1);
        CODE_TYPE nucleotide = categorial_sample<CODE_TYPE, ALPHABET_SIZE>(stationary);

        /* check for an outgroup */
        if (tree.outgroup()) {
                pt_generate_observations<CODE_TYPE, ALPHABET_SIZE>(
                        static_cast<const pt_leaf_t&>(*tree.outgroup()), stationary, nucleotide, observations);
        }
        /* traverse the treee */
        if (!tree.leaf()) {
                pt_generate_observations<CODE_TYPE, ALPHABET_SIZE>(
                        tree.left (), stationary, nucleotide, observations);
                pt_generate_observations<CODE_TYPE, ALPHABET_SIZE>(
                        tree.right(), stationary, nucleotide, observations);
        }
        return observations;
}

#endif /* PHYLOTREE_GENERATE_OBSERVATIONS_HH */
