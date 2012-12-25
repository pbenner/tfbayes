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
        yyparse();
        fclose(yyin);

        pt_root_t* pt_root = (pt_root_t*)pt_parsetree->convert();
        pt_parsetree->destroy();

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

const char* pt_leaf_name_rec(pt_node_t* node, size_t* count)
{
        const char* result = NULL;

        if (node->leaf()) {
                if (*count == 0) {
                        result = node->name.c_str();
                }
                else {
                        (*count)--;
                }
        }
        else {
                result = pt_leaf_name_rec(node->left, count);
                if (result == NULL) {
                        result= pt_leaf_name_rec(node->right, count);
                }
        }
        return result;
}

const char* pt_leaf_name(pt_node_t* node, size_t leaf)
{
        size_t count = leaf;

        return pt_leaf_name_rec(node, &count);
}

size_t pt_num_leafs(pt_node_t* node)
{
        size_t count = 0;

        if (node->leaf()) {
                count++;
        }
        else {
                count += pt_num_leafs(node->left);
                count += pt_num_leafs(node->right);
        }
        return count;
}

size_t pt_init(vector_t* observations, pt_node_t* node, size_t count)
{
        if (node->leaf()) {
                node->x = observations->vec[count];
                count++;
        }
        else {
                count = pt_init(observations, node->left,  count);
                count = pt_init(observations, node->right, count);
        }
        return count;
}

ssize_t pt_index_rec(pt_node_t* node, const string& name, size_t* count)
{
        ssize_t result = -1;

        if (node->leaf()) {
                if (node->name == name) {
                        result = *count;
                }
                (*count)++;
        }
        else {
                result = pt_index_rec(node->left,  name, count);
                if (result == -1) {
                        result = pt_index_rec(node->right, name, count);
                }
        }
        return result;
}

ssize_t pt_index(pt_root_t* pt_root, const char* name)
{
        size_t count = 0;
        string str(name);

        return pt_index_rec(pt_root, str, &count);
}

vector_t* pt_expectation(pt_root_t* pt_root, vector_t* observations, vector_t* prior)
{
        pt_init(observations, pt_root, 0);

        vector_t* result = alloc_vector(alphabet_size);
        pt_polynomial_t<code_t, alphabet_size> poly(pt_root);

        exponent_t<code_t, alphabet_size> alpha;
        for (size_t i = 0; i < alphabet_size; i++) {
                alpha[i] = prior->vec[i];
        }

        boost::array<double, alphabet_size> expectation =
//                pt_posterior_expectation_prime<code_t, alphabet_size>(poly, alpha);
                pt_posterior_expectation<code_t, alphabet_size>(poly, alpha);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = expectation[i];
        }

        return result;
}

vector_t* pt_approximate(pt_root_t* pt_root, vector_t* observations)
{
        pt_init(observations, pt_root, 0);

        vector_t* result = alloc_vector(alphabet_size);
        pt_polynomial_t<code_t, alphabet_size> poly(pt_root);

        pt_polynomial_t<code_t, alphabet_size> variational
                = dkl_approximate<code_t, alphabet_size>(poly);

        for (size_t i = 0; i < alphabet_size; i++) {
                result->vec[i] = variational.begin()->exponent()[i];
        }

        return result;
}

vector_t* pt_dkl_optimize(pt_root_t* pt_root, vector_t* observations)
{
        pt_init(observations, pt_root, 0);

        exponent_t<code_t, alphabet_size> alpha;
        alpha[0] = 1;
        alpha[1] = 1;
        alpha[2] = 1;
        alpha[3] = 1;
        alpha[4] = 1;

        vector_t* result = alloc_vector(alphabet_size);
        pt_polynomial_t<code_t, alphabet_size> poly(pt_root);

        pt_polynomial_t<code_t, alphabet_size> variational
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
