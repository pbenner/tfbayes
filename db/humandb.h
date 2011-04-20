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

#define HUMANDB_PAGE_SIZE 1024
#define HUMANDB_NUCLEOTIDES_PER_PAGE 128

void hdb_set_program_name(const char* program_name);
void hdb_free_program_name();
void hdb_set_error_file_pointer(FILE* error_file_pointer);
int  hdb_open(DB** _dbp, const char* db_name);
int  hdb_close(DB* dbp);
int  hdb_load_maf(DB *dbp, const char* maf);
int  hdb_get_sequence(DB *dbp, long pos_from, long pos_to, char* buf);

#endif /* HUMANDB_H */
