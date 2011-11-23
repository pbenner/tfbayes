/**
 *      @file  abstract_set.h
 *     @brief  Machine representation of the partitions of a finite set
 *
 *			   Partitions generator for a generic set of given size.
 *
 *    @author  P.-Y. Bourguignon (pyb), bourguig@mis.mpg.de
 *
 *  @internal
 *    Created  03/09/2010
 *   Revision  ---
 *    Company  Max-Planck Institute for Mathematics in the Sciences, Leipzig
 *  Copyright  Copyright (c) 2010, P.-Y. Bourguignon
 *
 * This source code is released for free distribution under the terms of the
 * GNU General Public License as published by the Free Software Foundation.
 */

#ifndef ABSTACT_SET_H
#define ABSTACT_SET_H

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

#define POW2(x) (1 << (x))
/**
 * @brief Represents an abstract symbol
 * 
 * Also used for the number of symbols (the alphabet size)
 */
typedef unsigned short symbol_t;

/**
 * @brief Represents an abstract subset
 * 
 * Bitwise decomposition of a subset_t indicates presence/absence of the symbols
 * This typing is valid for alphabets of size up to 8.
 */
typedef unsigned short subset_t;

/**
 * @brief Represents an abstract partition
 * 
 * The array of subset_t values is terminated by a '0' entry.
 */
typedef subset_t *partition_t;

/*
 * ! \brief Abstraction of the alphabet 
 */
typedef struct
{
  symbol_t    size;
  subset_t    nb_subsets;
  int         nb_partitions;
  partition_t *partitions;
  subset_t	  **containers;
} abstract_set_t;

/*
 * ! \brief Creation utility function 
 */
abstract_set_t *as_create (int alphabet_size);

/*
 * ! \brief Init utility function 
 */
int         as_init (abstract_set_t * as, int alphabet_size);

void        as_free (abstract_set_t * as);

int         as_subset_contains (subset_t subset, symbol_t symbol);

int         as_subset_size (subset_t subset);

subset_t    as_complement (const abstract_set_t * as, subset_t subset);

subset_t    as_append_symbol (subset_t subset, symbol_t symbol);

subset_t *	as_get_containers( const abstract_set_t * as, symbol_t symbol);

int 		as_get_largest_item(symbol_t alphabet_size, subset_t subset);

/*
 * ! \brief Return the number of components in a partition 
 */
int         as_partition_size (const partition_t partition);

void        as_append_subset (const abstract_set_t * as,
			      partition_t partition, subset_t subset);

int         _as_generate_partitions ();

__END_DECLS

#endif /* ABSTACT_SET_H */
