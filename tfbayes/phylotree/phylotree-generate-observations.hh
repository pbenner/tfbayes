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

#ifndef __TFBAYES_PHYLOTREE_PHYLOTREE_GENERATE_OBSERVATIONS_HH__
#define __TFBAYES_PHYLOTREE_PHYLOTREE_GENERATE_OBSERVATIONS_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cstdlib>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/statistics.hh>

/* AS: ALPHABET SIZE
 * AC: ALPHABET CODE TYPE
 * PC: POLYNOMIAL CODE TYPE
 */

template <size_t AS, typename AC>
void pt_generate_observations(const pt_node_t& node,
                              const std::vector<double>& stationary,
                              const AC parent_nucleotide,
                              std::vector<AC>& observations,
                              boost::random::mt19937& gen)
{
        double r = (double)rand()/RAND_MAX;
        AC nucleotide;
        if (r <= node.mutation_probability()) {
                nucleotide = categorical_sample<AS, AC>(stationary, gen);
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
                pt_generate_observations<AS, AC>(
                        node.left (), stationary, nucleotide, observations, gen);
                pt_generate_observations<AS, AC>(
                        node.right(), stationary, nucleotide, observations, gen);
        }
}

template <size_t AS, typename AC>
std::vector<AC>
pt_generate_observations(const pt_root_t& tree, const std::vector<double>& stationary,
                         boost::random::mt19937& gen)
{
        std::vector<AC> observations(tree.n_leaves, -1);
        AC nucleotide = categorical_sample<AS, AC>(stationary, gen);

        /* check for an outgroup */
        if (tree.outgroup()) {
                pt_generate_observations<AS, AC>(
                        static_cast<const pt_leaf_t&>(*tree.outgroup()), stationary, nucleotide, observations, gen);
        }
        /* traverse the treee */
        if (!tree.leaf()) {
                pt_generate_observations<AS, AC>(
                        tree.left (), stationary, nucleotide, observations, gen);
                pt_generate_observations<AS, AC>(
                        tree.right(), stationary, nucleotide, observations, gen);
        }
        return observations;
}

#endif /* __TFBAYES_PHYLOTREE_PHYLOTREE_GENERATE_OBSERVATIONS_HH__ */
