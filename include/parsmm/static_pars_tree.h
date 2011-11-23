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

#define GET_CHILD_OFFSET(as, node, subset) ((node)*(as)->nb_subsets+(subset))

#define GET_PARENT_OFFSET(as, node) (((node) - 1)/(as)->nb_subsets)

#define GET_SUBSET(as, node) ((subset_t)(((node) - 1) % (as)->nb_subsets) + 1)

#define GET_COUNTS(tree, node) ((tree)->counts + (tree)->as->size * (node))

#define GET_DIR_PARAMS(tree, node) ((tree)->dirichlet_params + (tree)->as->size * (node))


typedef struct static_pars_tree static_pars_tree_t ;

typedef unsigned long count_t ;

/*! Index of a node */
typedef unsigned long node_t ;

/*! Index of a context word (written in reverse) */
typedef unsigned long word_t ;

typedef enum {MV_UP, MV_DOWN, MV_SIBLING } move_t;

typedef enum {PRIOR_BD, PRIOR_CONSTANT} prior_family_t;

struct static_pars_tree {
	short depth ;
	const abstract_set_t * as ;
	node_t size ;
	count_t * counts ;
	double * scores ;
	double * dirichlet_params ;
};

static_pars_tree_t * pt_create(const abstract_set_t * as,
                               short depth);

void pt_free(static_pars_tree_t* tree);

unsigned long _pt_size(const abstract_set_t * as,
                       int depth);


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
double pt_ln_marginal_likelihood(const static_pars_tree_t * tree, 
                                 const count_t * obs);

double pt_ln_marginal_likelihood_add(const static_pars_tree_t * tree,
                                     const symbol_t ** contexts);

double pt_ln_marginal_likelihood_rm(const static_pars_tree_t * tree,
                                    const symbol_t ** contexts);

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

double _get_ln_score(const static_pars_tree_t * tree,
                     node_t node) ;

double _dirichlet_eval(const count_t * obs, 
                       const double* prior,
                       short alphabet_size) ;

void pt_set_prior(static_pars_tree_t * tree,
                  double weight,
                  prior_family_t prior_family);

int ipow(int base, int exp);

void sum_counts(count_t*, count_t*, count_t*, int); 

__END_DECLS

#endif /* STATIC_PARS_TREE_H */
