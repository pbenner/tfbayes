/* Copyright (C) 2010 Philipp Benner
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

#ifndef LINALG_H
#define LINALG_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stddef.h>

#include <gsl/gsl_matrix.h>

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

typedef struct {
        size_t size;
        double *vec;
} vector_t;

typedef struct {
        size_t rows;
        size_t columns;
        double **mat;
} matrix_t;

static inline
vector_t * alloc_vector(size_t size) {
        vector_t *v = (vector_t *)malloc(sizeof(vector_t));
        v->vec      = (double   *)malloc(sizeof(double) * size);
        v->size     = size;
        return v;
}

static inline
matrix_t * alloc_matrix(size_t rows, size_t columns) {
        matrix_t *m = (matrix_t *)malloc(sizeof(matrix_t));
        m->mat      = (double  **)malloc(sizeof(double *) * rows);
        m->rows     = rows;
        m->columns  = columns;
        size_t i;
        for (i = 0; i < rows; i++) {
                m->mat[i] = (double *)malloc(sizeof(double) * columns);
        }
        return m;
}

static inline
void free_vector(vector_t *v) {
        free(v->vec);
        free(v);
}

static inline
void free_matrix(matrix_t *m) {
        size_t i;
        for (i = 0; i < m->rows; i++) {
                free(m->mat[i]);
        }
        free(m->mat);
        free(m);
}

static inline
gsl_vector * to_gsl_vector(const vector_t *vector)
{
        gsl_vector *v = gsl_vector_alloc(vector->size);
        size_t i;

        for(i = 0; i < vector->size; i++) {
                gsl_vector_set(v, i, vector->vec[i]);
        }
        return v;
}

static inline
vector_t * from_gsl_vector(const gsl_vector * vector)
{
        size_t i;
        size_t size = vector->size;
        vector_t *v = alloc_vector(size);

        for(i = 0; i < size; i++) {
                v->vec[i] = gsl_vector_get(vector, i);
        }
        return v;
}

static inline
gsl_matrix * to_gsl_matrix(const matrix_t *matrix)
{
        gsl_matrix *m = gsl_matrix_alloc(matrix->rows, matrix->columns);
        size_t i, j;

        for(i = 0; i < matrix->rows; i++) {
                for(j = 0; j < matrix->columns; j++) {
                        gsl_matrix_set(m, i, j, matrix->mat[i][j]);
                }
        }
        return m;
}

static inline
matrix_t * from_gsl_matrix(const gsl_matrix * matrix)
{
        size_t i, j;
        size_t rows    = matrix->size1;
        size_t columns = matrix->size2;
        matrix_t *m    = alloc_matrix(rows, columns);

        for(i = 0; i < rows; i++) {
                for(j = 0; j < columns; j++) {
                        m->mat[i][j] = gsl_matrix_get(matrix, i, j);
                }
        }
        return m;
}

__END_DECLS

#endif /* LINALG_H */
