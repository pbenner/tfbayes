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

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>

#include <tfbayes/interface/common.hh>
#include <tfbayes/phylotree/interface.hh>
#include <tfbayes/phylotree/utility.hh>
#include <tfbayes/phylotree/phylotree-parser.h>

using namespace std;

__BEGIN_DECLS

// library interface
////////////////////////////////////////////////////////////////////////////////

extern FILE *yyin;

pt_root_t* pt_parse_file(const char* file_name)
{
        yyin = fopen(file_name, "r");
        if (yyin == NULL) {
                std_err(PERR, "Could not open phylogenetic tree");
        }

        list<pt_root_t*> tree_list = parse_tree_list(yyin);
        assert(tree_list.size() == 1);
        pt_root_t* pt_root = tree_list.front();
        fclose(yyin);

        return pt_root;
}

pt_node_t* pt_clone(pt_node_t* pt_node)
{
        return pt_node->clone();
}

void pt_destroy(pt_node_t* pt_node)
{
        pt_node->destroy();
}

char* pt_print(pt_root_t* pt_root)
{
        std::ostringstream stream;
        stream << pt_root;
        std::string str = stream.str();
        const char* chr = str.c_str();

        char* result = (char*)malloc((strlen(chr)+1)*sizeof(char));
        strcpy(result, chr);

        return result;
}

const char* pt_leaf_name(pt_root_t* root, size_t leaf)
{
        return (*root)(leaf)->name.c_str();
}

size_t pt_num_leafs(pt_root_t* root)
{
        return root->n_leafs;
}

ssize_t pt_index(pt_root_t* root, const char* name)
{
        return (*root)(name)->id;
}

vector_t* pt_expectation(pt_root_t* pt_root, vector_t* observations, vector_t* prior)
{
        vector_t* result = alloc_vector(alphabet_size);
        // convert observations to an std array
        vector<code_t> tmp(observations->vec, observations->vec+observations->size);
        // compute the polynomial
        polynomial_t<code_t, alphabet_size> poly = pt_polynomial<code_t, alphabet_size>(pt_root, tmp);

        exponent_t<code_t, alphabet_size> alpha;
        for (size_t i = 0; i < alphabet_size; i++) {
                alpha[i] = prior->vec[i];
        }

        boost::array<double, alphabet_size> expectation =
                pt_posterior_expectation<code_t, alphabet_size>(poly, alpha);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = expectation[i];
        }

        return result;
}

vector_t* pt_approximate(pt_root_t* pt_root, vector_t* observations)
{
        vector_t* result = alloc_vector(alphabet_size);
        // convert observations to an std array
        vector<code_t> tmp(observations->vec, observations->vec+observations->size);
        // compute the polynomial
        polynomial_t<code_t, alphabet_size> poly = pt_polynomial<code_t, alphabet_size>(pt_root, tmp);

        polynomial_t<code_t, alphabet_size> variational
                = dkl_approximate<code_t, alphabet_size>(poly);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = variational.begin()->exponent()[i];
        }

        return result;
}

vector_t* pt_dkl_optimize(pt_root_t* pt_root, vector_t* observations)
{
        vector_t* result = alloc_vector(alphabet_size);
        // convert observations to an std array
        vector<code_t> tmp(observations->vec, observations->vec+observations->size);
        // compute the polynomial
        polynomial_t<code_t, alphabet_size> poly = pt_polynomial<code_t, alphabet_size>(pt_root, tmp);

        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        alpha[4] = 1;

        polynomial_t<code_t, alphabet_size> variational
                = dkl_optimize(poly, alpha);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = variational.begin()->exponent()[i];
        }

        return result;
}

void pt_free(pt_root_t* pt_root)
{
        pt_root->destroy();
}

__END_DECLS
