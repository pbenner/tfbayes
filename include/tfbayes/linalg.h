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

#include <gsl/gsl_matrix.h>

typedef struct vector {
        int size;
        double *vec;
} Vector;

typedef struct matrix {
        int rows;
        int columns;
        double **mat;
} Matrix;

static inline
Vector * allocVector(int size) {
        Vector *v  = (Vector *)malloc(sizeof(Vector));
        v->vec     = (double *)malloc(sizeof(double) * size);
        v->size    = size;
        return v;
}

static inline
Matrix * allocMatrix(int rows, int columns) {
        Matrix *m  = (Matrix * )malloc(sizeof(Matrix));
        m->mat     = (double **)malloc(sizeof(double *) * rows);
        m->rows    = rows;
        m->columns = columns;
        int i;
        for (i = 0; i < rows; i++) {
                m->mat[i] = (double *)malloc(sizeof(double) * columns);
        }
        return m;
}

static inline
void freeVector(Vector *v) {
        free(v->vec);
        free(v);
}

static inline
void freeMatrix(Matrix *m) {
        int i;
        for (i = 0; i < m->rows; i++) {
                free(m->mat[i]);
        }
        free(m->mat);
        free(m);
}

static inline
gsl_vector * toGslVector(const Vector *vector)
{
        gsl_vector *v = gsl_vector_alloc(vector->size);
        int i;

        for(i = 0; i < vector->size; i++) {
                gsl_vector_set(v, i, vector->vec[i]);
        }
        return v;
}

static inline
Vector * fromGslVector(const gsl_vector * vector)
{
        int i;
        int size    = vector->size;
        Vector *v   = allocVector(size);

        for(i = 0; i < size; i++) {
                v->vec[i] = gsl_vector_get(vector, i);
        }
        return v;
}

static inline
gsl_matrix * toGslMatrix(const Matrix *matrix)
{
        gsl_matrix *m = gsl_matrix_alloc(matrix->rows, matrix->columns);
        int i, j;

        for(i = 0; i < matrix->rows; i++) {
                for(j = 0; j < matrix->columns; j++) {
                        gsl_matrix_set(m, i, j, matrix->mat[i][j]);
                }
        }
        return m;
}

static inline
Matrix * fromGslMatrix(const gsl_matrix * matrix)
{
        int i, j;
        int rows    = matrix->size1;
        int columns = matrix->size2;
        Matrix *m   = allocMatrix(rows, columns);

        for(i = 0; i < rows; i++) {
                for(j = 0; j < columns; j++) {
                        m->mat[i][j] = gsl_matrix_get(matrix, i, j);
                }
        }
        return m;
}

#endif /* LINALG_H */
