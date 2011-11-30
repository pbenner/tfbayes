/*! \file static_pars_tree.h
 *  \brief Parsimonious Context Tree (static version)
 *
 * <H3>Indexing rationale</H3>
 *  <UL>
 *   <LI>SYMBOLS are numbered from 0 to alphabet_size-1</LI>
 *
 *   <LI>SUBSETS are numbered from 0 to nb_subsets-1, 
   with bit i set to 1 <=> symbol i in subset.</LI>
 * 
 *   <LI>NODES in the tree are identified by an integer n,
 * arrays are looked up using n * CELL_SIZE.</LI>
 *  </UL>
 */

#ifndef STATIC_PARS_TREE_H
#define STATIC_PARS_TREE_H

#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#include <parsmm/abstract_set.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct static_pars_tree static_pars_tree_t;

typedef unsigned long count_t;
typedef unsigned short depth_t;

/*! Index of a node */
typedef unsigned long node_t;

/*! Index of a context word (written in reverse) */
typedef unsigned long word_t;

typedef enum {MV_UP, MV_DOWN, MV_SIBLING} move_t;
typedef enum {PRIOR_BD, PRIOR_CONSTANT} prior_family_t;

struct static_pars_tree {
	depth_t depth;
	const abstract_set_t * as;
	node_t size;
	count_t * counts;
	double * scores;
	double * dirichlet_params;
        /* The array of pointers over nodes in the pars_tree */
        node_t * node_ids;
        node_t * new_node_ids;
};

static_pars_tree_t * pt_create(const abstract_set_t * as,
                               depth_t depth);

void pt_init(static_pars_tree_t* pt);
void pt_free(static_pars_tree_t* tree);


/*! \brief Compute the marginal likelihood for counts and prior provided
 *  \sa pt_get_prior()
 *
 *  \warning Counts are stored in the tree, subsequent calls on the same 
 *  tree object overwrite them.
 *
 *  Further updates of the count statistics
 *  and recomputations of the marginal likelihood should use the variants
 *  suffixed \c add and \c rm.
 *  	
 *  \arg \c prior An array of size tree->size * alphabet_size, providing the 
 * 		 Dirichlet parameters for each node.
 *
 *  \arg \c obs An array of size alphabet_size^(tree->depth+1) providing 
 * 		 the count statistics (1 cell with the alphabet_size counts stored
 *  	 consecutively, cells sorted by alphabetical order of the context
 *		 read from right to left).
 *
 *  \return the logarithm of the marginal likelihood.
 */
double pt_ln_marginal_likelihood(static_pars_tree_t * tree, 
                                 const count_t * obs);

/*! \brief Prepares a prior array for use in pt_marginal_likelihood()
 *  <UL>
 *  <LI> \c prior_type Choose among:
 *  	<UL>
 *		<LI> \c PRIOR_BD the BDeu prior of Heckermann & co.</LI>
 *		<LI> \c PRIOR_CONSTANT the constant prior</LI>
 *		</UL>
 *	</LI>
 *  </UL>
 */

void pt_set_prior(static_pars_tree_t * tree,
                  double weight,
                  prior_family_t prior_family);

__END_DECLS

#endif /* STATIC_PARS_TREE_H */
