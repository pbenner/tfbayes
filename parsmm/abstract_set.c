/*
 *      @file  abstract_set.c
 *      @brief  Helper functions manipulating abstract_sets
 *
 *	Definitions of the functions declared in abstract_set.h
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <parsmm/abstract_set.h>

abstract_set_t * 
as_create(int alphabet_size)
{
	abstract_set_t *as ;
	if ((as = (abstract_set_t*)malloc(sizeof(abstract_set_t))))
	{
		if ( as_init(as, alphabet_size) )
			return as ;
		else
			return NULL ;
	}
	else
		return NULL ;
}

int
as_init(
        abstract_set_t *as, 
        int alphabet_size )
{
	symbol_t ii ;
	subset_t jj ;
	subset_t umask ;
	as->size = alphabet_size ;
	as->nb_subsets = POW2(alphabet_size) -1 ;
	as->containers = (subset_t **) malloc( sizeof(subset_t *)*alphabet_size ) ;
	for ( ii = 0 ; ii < alphabet_size ; ii++ )
	{
		as->containers[ii] = (subset_t *) malloc( sizeof(subset_t) * (1 << ( alphabet_size - 1 ))) ;
		umask = ( (subset_t) ~0 ) << ii ;
		for ( jj = 0 ; jj < (as->nb_subsets + 1)/2; jj++ )
		{
			as->containers[ii][jj] = (( jj & umask ) << 1 ) | ( jj & ~umask ) | (1 << ii) ;
		}
	}
	if (_as_generate_partitions(as))
		return 1 ;
	else
		return 0 ;
}

void
as_free(abstract_set_t * as)
{
	int ii = 0 ;
	while (as->partitions[ii])
	{
		free(as->partitions[ii]) ;
		ii++ ;
	}
	free(as->partitions) ;
	for (ii = 0 ; ii < as->size ; ii++ )
	{
		free(as->containers[ii]) ;
	}
	free(as->containers) ;
	free(as) ;
}

int 
as_subset_contains(subset_t subset, symbol_t symbol)
{
	return ( subset >> symbol ) & 1 ;
}

subset_t 
as_append_symbol(subset_t subset, symbol_t symbol)
{
	return subset | (1 << symbol) ;
}

int 
as_subset_size(subset_t subset)
{
	return __builtin_popcount(subset) ;
}

subset_t 
as_subset_complement(
        const abstract_set_t *as, 
        subset_t subset)
{
	return (1 << as->size) - 1 - subset ;
}

/* Partitions */


int
as_partition_size(
        const partition_t partition)
{
	int num = 0 ;
	while (partition[num])
	{
		num += 1 ;
	}
	return num ;	
}

void
as_append_subset(
        const abstract_set_t *as, 
        partition_t partition, 
        subset_t subset )
{
	int last = 0 ;
	while (partition[last] && ( last < as->size - 1 ))
	{
		last += 1 ;
	}
	partition[last] = subset ;
	partition[last+1] = 0 ;
}

int 
_as_generate_partitions(abstract_set_t *as)
{
	partition_t * previous = NULL ;
	partition_t * current  = NULL ;
	int alphabet_size = 1 ;
	int next_size = 0 ;
	int ii ;
	int jj = 0 ;
	int kk ;
/*  	subset_t subset ;*/
	#ifdef DEBUG
		printf("*** Starting partitions generation ***\n") ;
	#endif /* DEBUG */

	/* Init (alphabet_size = 1) */
	current = malloc( 2 * sizeof( partition_t ) ) ;
	current[0] = malloc( 2 * sizeof( subset_t ) ) ;
	current[0][0] = 1 ;
	current[0][1] = 0 ;
	current[1] = NULL ;
	next_size = 2 ;

	/* Iterations */
	
	for ( alphabet_size = 2 ; alphabet_size < as->size + 1 ; alphabet_size++ )
	{
		#ifdef DEBUG
		printf("In generate_partitions: Alphabet size: %d\n", alphabet_size) ;
		#endif /* DEBUG */
		/* Rotate the arrays */
		
		previous = current ;
		current = malloc( ( next_size + 1 ) * sizeof( partition_t ) ) ;
		
		/* Reset utility variables */
		ii = 0 ; /* used to iterate over the previous partitions */
		jj = 0 ; /* used to index current partitions */
		next_size = 0 ; 
		/* Iterate on the previous partitions */
		while ( previous[ii] )
		{
			kk = 0 ; /* used to iterate over the subsets in the previous partition */
			
			/* Add the new letter as a singleton */
			current[jj] = malloc( (alphabet_size + 1) * sizeof( subset_t ) ) ;
			memcpy( current[jj], previous[ii], alphabet_size * sizeof( subset_t ) ) ;
			as_append_subset( as, current[jj], POW2( alphabet_size - 1 ) ) ;
			next_size += 1 + as_partition_size(current[jj]) ;
			jj += 1 ;
			/* Create one partition by adding the new element to each subset */
			while ( previous[ii][kk] )
			{
				current[jj] = malloc( (alphabet_size + 1) * sizeof( subset_t ) ) ;
				memcpy( current[jj], previous[ii], alphabet_size * sizeof( subset_t ) ) ;
				current[jj][kk] = as_append_symbol( current[jj][kk], alphabet_size - 1) ;
				next_size += 1 + as_partition_size(current[jj]) ;
				jj += 1 ;
				kk += 1 ;
			}
			ii += 1 ;
		}
		current[jj] = NULL ;
		ii = 0 ;
		while ( previous[ii] )
		{
			free( previous[ii] ) ;
			ii += 1 ;
		}
		free(previous) ;	 
	}
	as->partitions = current ;	
	as->nb_partitions = jj ;

	#ifdef DEBUG
		printf("*** Partitions generation done ***\n") ;
	#endif /* DEBUG */

	return as->nb_partitions ;
}

#if defined DEBUG || defined TEST
int partition_check(
        const abstract_set_t *as, 
        const partition_t partition)
{
	subset_t sum = 0 ;
	int ii = 0 ;
	while( partition[ii] )
	{
		sum = sum | partition[ii] ;
		ii += 1 ;
		if ( ( sum & partition[ii] ) != 0 )
		{
			printf("Overlap between subsets: sum is %u, subset is %u\n", sum, partition[ii]) ;
			return 0 ;
		}
	}
	if ( sum == (subset_t) ( (1 << as->size) - 1 ) )
	{
		return 1 ;
	}
	else
	{
		printf("Subsets do not cover the alphabet\n") ;
		return 0 ;
	}
}
#endif /* DEBUG || TEST */

int as_get_largest_item(symbol_t alphabet_size, subset_t subset)
{
	/* WARNING: Returns 0 also for an empty subset */
	int ii = alphabet_size - 1; 
	do {
		if ( ((1 << ii) & subset ) == (1 << ii))
			break ;
	}
	while (--ii) ;
	return ii ;
}
