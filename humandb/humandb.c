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

#include "humandb.h"

static char* _program_name;
static FILE* _error_file_pointer;

static u_int32_t num_records(DB* dbp)
{
        u_int32_t num;
        DB_BTREE_STAT* sp = NULL;
        dbp->stat(dbp, NULL, &sp, DB_FAST_STAT);
        num = sp->bt_ndata;

        free(sp);

        return num;
}

void hdb_set_program_name(const char* program_name)
{
        _program_name = (char *)malloc(sizeof(char)*strlen(program_name)+1);
        strcpy(_program_name, program_name);
}

void hdb_free_program_name()
{
        free(_program_name);
}

void hdb_set_error_file_pointer(FILE* error_file_pointer)
{
        _error_file_pointer = error_file_pointer;
}

int hdb_open(DB** _dbp, const char* db_file_name, const char* db_name)
{
        DB* dbp;
        int ret;
        u_int32_t open_flags;

        ret = db_create(&dbp, NULL, 0);
        if (ret != 0) {
                fprintf(stderr, "%s: %s\n", _program_name,
                        db_strerror(ret));
                return ret;
        }

        /* Set up error handling for this database */
        dbp->set_errfile(dbp, _error_file_pointer);
        dbp->set_errpfx(dbp, _program_name);

        open_flags = DB_CREATE;

        /* page size in bytes */
        dbp->set_pagesize(dbp, HUMANDB_PAGE_SIZE);
        /* record length in bytes */
        dbp->set_re_len(dbp, HUMANDB_RECORD_LENGTH);
        dbp->set_re_pad(dbp, '\0');
        ret = dbp->open(dbp,          /* Pointer to the database */
                        NULL,         /* Txn pointer */
                        db_file_name, /* File name */
                        db_name,      /* Logical db name */
                        DB_RECNO,     /* Database type (using recno) */
                        open_flags,   /* Open flags */
                        0);           /* File mode. Using defaults */

        if (ret != 0) {
                dbp->err(dbp, ret, "Database '%s:%s' open failed.", db_file_name, db_name);
                return ret;
        }

        (*_dbp) = dbp;
        return 0;
}

int hdb_close(DB* dbp)
{
        int ret;

        ret = dbp->sync(dbp, 0);
        ret = dbp->close(dbp, 0);

        if (ret != 0) {
                fprintf(_error_file_pointer, "Database close failed: %s\n",
                        db_strerror(ret));
        }
        return ret;
}

int hdb_load_maf(DB *dbp, const char* maf)
{
        DBT key, data;
        char buf[HUMANDB_RECORD_LENGTH+1];
        int ret;

        if (num_records(dbp) > 0) {
                fprintf(_error_file_pointer, "Database has already been initialized\n");
                return -1;
        }

        FILE* ifp = fopen(maf, "r");
        if (ifp == NULL) {
                fprintf(_error_file_pointer, "Error opening file '%s'\n", maf);
                return -1;
        }

        db_recno_t page = 1;
        while (fgets(buf, HUMANDB_RECORD_LENGTH+1, ifp) != NULL) {
                size_t len = strlen(buf);

                if (buf[len-1] == '\n') {
                        len--;
                }

                memset(&key,  0, sizeof(DBT));
                memset(&data, 0, sizeof(DBT));

                key.data  = &page;
                key.size  = sizeof(db_recno_t);

                data.data = buf;
                data.size = len;

                ret = dbp->put(dbp, 0, &key, &data, 0);
                if (ret != 0) {
                        dbp->err(dbp, ret, "Error inserting: '%s'", buf);
                        return ret;
                }
                page++;
        }

        dbp->sync(dbp, 0);
        fclose(ifp);
        return 0;
}

int hdb_get_sequence(DB *dbp, long pos_from, long n_nucleotides, char* buf)
{
        DBT key, data;
        long pos;
        db_recno_t from_page, from_page_offset; 
        db_recno_t i_page, i_page_offset;

        memset(&key,  0, sizeof(DBT));
        memset(&data, 0, sizeof(DBT));

        /* logical record numbers start at 1 */
        from_page        = pos_from/HUMANDB_RECORD_LENGTH + 1;
        from_page_offset = pos_from%HUMANDB_RECORD_LENGTH;

        /* loop through the database */
        pos           = 0;
        i_page_offset = from_page_offset;
        for(i_page = from_page; pos < n_nucleotides; i_page++) {
                /* initialize key */
                key.data = &i_page;
                key.size = sizeof(db_recno_t);

                /* retreive the next page */
                if (dbp->get(dbp, NULL, &key, &data, 0) != 0) {
                        goto err;
                }

                while (pos < n_nucleotides && i_page_offset < data.size)
                {
                        buf[pos++] = ((char *)data.data)[i_page_offset++];
                }
                /* reset pointer */
                i_page_offset = 0;
        }
err:
        buf[pos] = '\0';
        return 0;
}

int hdb_get_sequence_pure(DB *dbp, long pos_from, long n_nucleotides, char* buf)
{
        DBT key, data;
        long pos;
        db_recno_t from_page, from_page_offset; 
        db_recno_t i_page, i_page_offset;

        memset(&key,  0, sizeof(DBT));
        memset(&data, 0, sizeof(DBT));

        /* logical record numbers start at 1 */
        from_page        = pos_from/HUMANDB_RECORD_LENGTH + 1;
        from_page_offset = pos_from%HUMANDB_RECORD_LENGTH;

        /* loop through the database */
        pos           = 0;
        i_page_offset = from_page_offset;
        for(i_page = from_page; pos < n_nucleotides; i_page++) {
                /* initialize key */
                key.data = &i_page;
                key.size = sizeof(db_recno_t);

                /* retreive the next page */
                if (dbp->get(dbp, NULL, &key, &data, 0) != 0) {
                        goto err;
                }

                while (pos < n_nucleotides && i_page_offset < data.size)
                {
                        if (is_nucleotide(((char *)data.data)[i_page_offset])) {
                                buf[pos++] = ((char *)data.data)[i_page_offset++];
                        }
                        else {
                                i_page_offset++;
                        }
                }
                /* reset pointer */
                i_page_offset = 0;
        }
err:
        buf[pos] = '\0';
        return 0;
}
