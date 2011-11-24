/*
 *      @file  static_pars_tree.c
 *      @brief  Parsimonious context tree
 *
 *     @author  P.-Y. Bourguignon (pyb), bourguig@mis.mpg.de
 *
 *   @internal
 *     Created  03/09/2010
 *     Revision  ---
 *     Company  Max-Planck Institute for Mathematics in the Sciences, Leipzig
 *     Copyright (c) 2010, P.-Y. Bourguignon
 *
 * This source code is released for free distribution under the terms of the
 * GNU General Public License as published by the Free Software Foundation.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <parsmm/static_pars_tree.h>
#include <tfbayes/exception.h>
#include <tfbayes/logarithmetic.h>

static
int ipow(int base, int exp) 
{
        int result = 1;
        while (exp)
        {
                if (exp & 1) {
                        result *= base;
                }
                exp >>= 1;
                base *= base;
        }
        return result;
} 

static
void sum_counts(count_t* dest, count_t* src1, count_t* src2, int size) 
{
        int ii;
        for (ii = 0; ii < size; ii++) {
                dest[ii] = src1[ii] + src2[ii];
        }
}

static
unsigned long pt_size(
        const abstract_set_t * as,
        int depth)
{
        int ii;
        unsigned long powk = as->nb_subsets;
        for (ii = 0; ii < depth; ii++)
        {
                powk = powk * as->nb_subsets;
        }
        return (powk - 1) / (as->nb_subsets - 1);
}

static_pars_tree_t * pt_create(
        const abstract_set_t * as,
        short depth)
{
        static_pars_tree_t * tree;
        tree = (static_pars_tree_t *) malloc(sizeof(static_pars_tree_t));
        tree->as = as;
        tree->size = pt_size(as, depth);
        tree->depth = depth;
        tree->counts = (count_t *)malloc(sizeof(unsigned long) * tree->size * as->size);
        tree->scores = (double *)malloc(sizeof(double) * tree->size);
        tree->dirichlet_params = (double *)malloc(sizeof(double) * tree->size * as->size);
        return tree;
}

void pt_free(static_pars_tree_t * tree)
{
        free(tree->dirichlet_params);
        free(tree->scores);
        free(tree->counts);
        free(tree);
}

static
double dirichlet_eval(
        const count_t * obs,
        const double * dirichlet_params,
        short alphabet_size)
{
        double sum_params = 0.0;
        double sum_counts = 0.0;
        double result = 0.0;
        int ii;

        print_debug("[dirichlet_eval]\n");

        for ( ii = 0; ii < alphabet_size; ii++ )
        {
                sum_params += dirichlet_params[ii];
                sum_counts += obs[ii];
                result += gsl_sf_lngamma(dirichlet_params[ii] + obs[ii]) - gsl_sf_lngamma(dirichlet_params[ii]);
                print_debug("%f,%lu\n", dirichlet_params[ii], obs[ii]);
        }
        result += gsl_sf_lngamma(sum_params) - gsl_sf_lngamma(sum_params + sum_counts);
        print_debug("Result: %e\n", result);

        return result;
}

static
double get_ln_score(
        const static_pars_tree_t * tree,
        node_t node)
{
        double result = 0;
        double partial_sum;
        short ii, jj; /* partition index, subset index */
        partition_t partition;

        print_debug("[get_ln_score node: %lu]\n",node);
        if ( node * tree->as->nb_subsets + 2 > tree->size )
        {
                /* Node is a leaf */
                result = dirichlet_eval(GET_COUNTS(tree, node), GET_DIR_PARAMS(tree, node), tree->as->size);
        }
        else
        {
                /* Node is internal */
                for ( ii = 0; ii < tree->as->nb_partitions; ii++ )
                {
                        partition = tree->as->partitions[ii];
                        
                        jj = 0;
                        partial_sum = 0;
                        
                        while ( partition[jj] != 0 )
                        {
                                partial_sum += tree->scores[GET_CHILD_OFFSET(tree->as, node, partition[jj])];
                                jj++; 
                        }
                        
                        result = ii == 0 ? partial_sum : logadd(result, partial_sum);
                }
        }
        return result;
}

double pt_ln_marginal_likelihood(
        const static_pars_tree_t* tree,
        const count_t * obs)
{
        /* The running index over the contexts */
        word_t context_id = 0;

        symbol_t symbol = 0; /* TOCHECK Maintain this upon iteration */

        /* The array of pointers over nodes in the pars_tree */
        node_t * node_ids ;
        node_t * new_node_ids;

        void * pivot;

        unsigned long ii, jj, kk;

        /* The number of nodes in the counts array */
        unsigned long context_id_max;

        /* The number of nodes containing a symbol */
        unsigned long container_max;

        /* Record the last move */
        move_t last_move;

        /* INITIALIZATION */
        context_id_max = ipow(tree->as->size, tree->depth);
        container_max = ( tree->as->nb_subsets + 1 ) / 2;
        /* First blank the counts and scores array */
        memset(tree->counts, 0, sizeof(unsigned long) * tree->size * tree->as->size);
        memset(tree->scores, 0, sizeof(double) * tree->size);
        /* Allocate memory for the stack of pointers 
         * Exactly half of the subsets contain a given symbol
         */
        node_ids = (node_t *) malloc(sizeof(node_t) * (1 << ((tree->as->size)*tree->depth)));
        new_node_ids = (node_t *) malloc(sizeof(node_t) * (1 << ((tree->as->size) * tree->depth)));
        memset(node_ids, (unsigned int) 0, sizeof(node_t) * (1 << ((tree->as->size)*tree->depth)));
        node_ids[0] = 0;
        last_move = MV_DOWN; /* Not to perturb the first iteration */

        /* ITERATIONS */
        do {
                memset(new_node_ids, (unsigned int) 0, sizeof(node_t) * (1 << ((tree->as->size)*tree->depth)));

                if (last_move != MV_UP)
                {
                        /* Copy and compute counts */
                        ii = 0;
                        do {
                                if (as_subset_size(GET_SUBSET(tree->as, node_ids[ii])) <= 1) {
                                        memcpy(GET_COUNTS(tree, node_ids[ii]), obs + context_id * tree->as->size, tree->as->size * sizeof(count_t));
                                        print_debug("Copying counts - src: %lu\tdst: %lu\n", context_id, node_ids[ii]);
                                }
                                else if ( GET_SUBSET(tree->as, node_ids[ii]) < (1 << (symbol + 1)))
                                {
                                        print_debug("Summing counts - current: %lu\n", node_ids[ii]);
                                        sum_counts(GET_COUNTS(tree, node_ids[ii]),GET_COUNTS(tree, node_ids[ii] - (1 << symbol)), GET_COUNTS(tree, GET_CHILD_OFFSET(tree->as,GET_PARENT_OFFSET(tree->as,node_ids[ii]),(1 << symbol))),tree->as->size);
                                }
                                if ( context_id * tree->as->size + 2 > context_id_max )
                                        tree->scores[node_ids[ii]] = get_ln_score(tree, node_ids[ii]);
                                ii++;
                        }
                        while (node_ids[ii]);

                                
                }
                else {
                        ii = 0;
                        do {
                                print_debug("[pt_ln_marginal_likelihood] Node: %lu - Computing score\n", node_ids[ii]);
                                tree->scores[node_ids[ii]] = get_ln_score(tree,node_ids[ii]);
                                ii++;
                        }
                        while ( node_ids[ii] != 0 );
                }

                /* Move the pointer over observations */
                if (( context_id * tree->as->size + 1 < context_id_max ) &&
                    ( last_move != MV_UP ))
                {
                        /* GO DOWN */
                        context_id = context_id * tree->as->size + 1;
                        symbol = 0;
                        last_move = MV_DOWN;
                        ii = 0;
                        kk = 0;
                        do {
                                for ( jj = 0; jj < container_max; jj++ )
                                        new_node_ids[kk++] = GET_CHILD_OFFSET(tree->as, node_ids[ii],tree->as->containers[0][jj]);
                                ii++;
                        }
                        while ( node_ids[ii] != 0 );
                }

                else if ( context_id % tree->as->size != 0 )
                {
                        /* GO TO NEXT SIBLING */
                        context_id = context_id + 1;
                        symbol++;
                        last_move = MV_SIBLING;
                        ii = 0;
                        kk = 0;
                        do {
                                for ( jj = 0; jj < container_max; jj++ )
                                        new_node_ids[kk++] = GET_CHILD_OFFSET(tree->as, GET_PARENT_OFFSET(tree->as,node_ids[ii]), tree->as->containers[(context_id - 1) % 4][jj]);
                                ii += container_max;
                        }
                        while ( node_ids[ii] != 0 );
                }

                else
                {
                        /* GO UP */
                        context_id = (context_id - 1)/tree->as->size;
                        symbol = (context_id - 1) % tree->as->size;
                        last_move = MV_UP;
                        ii = 0;
                        kk = 0;
                        do {
                                new_node_ids[kk++] = GET_PARENT_OFFSET(tree->as, node_ids[ii]);
                                ii += container_max;
                        }
                        while ( node_ids[ii] != 0 );
                }
                pivot = node_ids;
                node_ids = new_node_ids;
                new_node_ids = pivot;
        }
        while ( context_id != 0 );

        tree->scores[0] = get_ln_score(tree,0);

        /* CLEANUP */
        free(node_ids);
        free(new_node_ids);

        return tree->scores[0];
}
