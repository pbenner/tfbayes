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

#ifndef HUMANDB_H
#define HUMANDB_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tfbayes/linalg.h>

#define HUMANDB_PAGE_SIZE 1024
#define HUMANDB_RECORD_LENGTH 128

void hdb_set_num_threads(size_t num_threads);
void hdb_set_program_name(const char* program_name);
void hdb_free_program_name();
void hdb_set_error_file_pointer(FILE* error_file_pointer);
int  hdb_open(DB** _dbp, const char* db_file_name, const char* db_name, u_int32_t open_flags);
int  hdb_close(DB* dbp);
int  hdb_load_maf(DB *dbp, const char* maf);
int  hdb_get_sequence(DB *dbp, long pos_from, long n_nucleotides, char* buf);
int  hdb_get_sequence_pure(DB *dbp, long pos_from, long n_nucleotides, char* buf);
int  hdb_search(DB* dba[], int dba_n, const char* db_names[], const char* pattern, int pattern_n);
int  hdb_search_pwm(DB* dbp_list[], const int dbp_list_n, const char* db_names[], const Matrix* pwm, const double threshold);

static inline
int is_nucleotide(char S) {
        if (S == 'A' || S == 'C' || S == 'G' || S == 'T' ||
            S == 'a' || S == 'c' || S == 'g' || S == 't') {
                return 1;
        }
        else {
                return 0;
        }
}

#endif /* HUMANDB_H */
