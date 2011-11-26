/* Copyright (C) 2011 Philipp Benner
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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <db.h>

#include <tfbayes/linalg.h>

#include "humandb.h"
#include "statistics-codons.h"
#include "interface.h"

// export functions for interface.py
vector_t * _alloc_vector(size_t size) {
        return alloc_vector(size);
}
matrix_t * _alloc_matrix(size_t rows, size_t columns) {
        return alloc_matrix(rows, columns);
}
void _free_vector(vector_t *v) { free_vector(v); }
void _free_matrix(matrix_t *m) { free_matrix(m); }
void _free(void *ptr)          { free(ptr); }

void _hdb_init(const char* program_name, int num_threads)
{
        hdb_set_num_threads((size_t)num_threads);
        hdb_set_program_name(program_name);
        hdb_set_error_file_pointer(stderr);
}

void _hdb_free()
{
        hdb_free_program_name();
}

DB* _hdb_create(const char* db_file_name, const char* db_name)
{
        DB* dbp;

        hdb_open(&dbp, db_file_name, db_name, DB_CREATE);

        return dbp;
}

DB* _hdb_open_ro(const char* db_file_name, const char* db_name)
{
        DB* dbp;

        hdb_open(&dbp, db_file_name, db_name, DB_RDONLY);

        return dbp;
}

DB* _hdb_open(const char* db_file_name, const char* db_name)
{
        DB* dbp;

        hdb_open(&dbp, db_file_name, db_name, 0);

        return dbp;
}

void _hdb_close(DB* dbp)
{
        hdb_close(dbp);
}

void _hdb_load_maf(DB* dbp, const char* maf)
{
        hdb_load_maf(dbp, maf);
}

void _hdb_get_sequence(DB* dbp, size_t pos, size_t n_nucleotides, char* buf) {
        hdb_get_sequence(dbp, pos, n_nucleotides, buf);
}

void _hdb_get_sequence_pure(DB* dbp, size_t pos, size_t n_nucleotides, char* buf) {
        hdb_get_sequence_pure(dbp, pos, n_nucleotides, buf);
}

void _hdb_search_pwm(DB* dbp_list[], const int dbp_list_n, const char* db_names[], const matrix_t* pwm, double threshold) {
        hdb_search_pwm(dbp_list, dbp_list_n, db_names, pwm, threshold);
}

void _hdb_search(DB* dbp_list[], size_t dbp_list_n, const char* db_names[], const char* sequence) {
        int sequence_n = strlen(sequence);
        hdb_search(dbp_list, dbp_list_n, db_names, sequence, sequence_n);
}

void _hdb_count_codons(DB *dbp, long positions[], size_t n, long result[]) {
        hdb_count_codons(dbp, positions, n, result);
}

void _hdb_count_codons_upstream(DB *dbp, size_t len, long result[]) {
        hdb_count_codons_upstream(dbp, len, result);
}
