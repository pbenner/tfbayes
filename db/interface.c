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
#include "interface.h"

// export functions for interface.py
Vector * _allocVector(int size)              { return allocVector(size); }
void     _freeVector(Vector *v)              { freeVector(v); }
Matrix * _allocMatrix(int rows, int columns) { return allocMatrix(rows, columns); }
void     _freeMatrix(Matrix *m)              { freeMatrix(m); }
void     _free(void *ptr)                    { free(ptr); }

void _hdb_init(const char* program_name)
{
        hdb_set_program_name(program_name);
        hdb_set_error_file_pointer(stderr);
}

void _hdb_free()
{
        hdb_free_program_name();
}

DB* _hdb_open(const char* db_file_name, const char* db_name)
{
        DB* dbp;

        hdb_open(&dbp, db_file_name, db_name);

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
        hdb_get_sequence(dbp, pos, pos+n_nucleotides-1, buf);
}
