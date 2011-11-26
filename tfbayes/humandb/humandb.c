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
#include <ctype.h>

#include <pthread.h>
#include <db.h>

#include <tfbayes/exception.h>

#include "humandb.h"

static size_t _num_threads;
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

void hdb_set_num_threads(size_t num_threads)
{
        _num_threads = num_threads;
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

// open & close
////////////////////////////////////////////////////////////////////////////////

int hdb_open(DB** _dbp, const char* db_file_name, const char* db_name, u_int32_t open_flags)
{
        DB* dbp;
        int ret;

        ret = db_create(&dbp, NULL, 0);
        if (ret != 0) {
                fprintf(stderr, "%s: %s\n", _program_name,
                        db_strerror(ret));
                return ret;
        }

        /* Set up error handling for this database */
        dbp->set_errfile(dbp, _error_file_pointer);
        dbp->set_errpfx(dbp, _program_name);

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
                (*_dbp) = NULL;
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

// obtain sequences
////////////////////////////////////////////////////////////////////////////////

int hdb_get_sequence(DB *dbp, long pos_from, long n_nucleotides, char* buf)
{
        DBT key, data;
        long pos;
        db_recno_t from_rec, from_rec_offset; 
        db_recno_t i_rec, i_rec_offset;
        int ret;

        memset(&key,  0, sizeof(DBT));
        memset(&data, 0, sizeof(DBT));

        /* logical record numbers start at 1 */
        from_rec        = pos_from/HUMANDB_RECORD_LENGTH + 1;
        from_rec_offset = pos_from%HUMANDB_RECORD_LENGTH;

        /* loop through the database */
        pos          = 0;
        i_rec_offset = from_rec_offset;
        for(i_rec = from_rec; pos < n_nucleotides; i_rec++) {
                /* initialize key */
                key.data = &i_rec;
                key.size = sizeof(db_recno_t);

                /* retreive the next rec */
                if ((ret = dbp->get(dbp, NULL, &key, &data, 0)) != 0) {
                        buf[pos] = '\0';
                        return ret;
                }

                while (pos < n_nucleotides && i_rec_offset < data.size)
                {
                        buf[pos++] = ((char *)data.data)[i_rec_offset++];
                }
                /* reset pointer */
                i_rec_offset = 0;
        }

        return 0;
}

int hdb_get_sequence_pure(DB *dbp, long pos_from, long n_nucleotides, char* buf)
{
        DBT key, data;
        long pos;
        db_recno_t from_rec, from_rec_offset; 
        db_recno_t i_rec, i_rec_offset;
        int ret;

        memset(&key,  0, sizeof(DBT));
        memset(&data, 0, sizeof(DBT));

        /* logical record numbers start at 1 */
        from_rec        = pos_from/HUMANDB_RECORD_LENGTH + 1;
        from_rec_offset = pos_from%HUMANDB_RECORD_LENGTH;

        /* loop through the database */
        pos          = 0;
        i_rec_offset = from_rec_offset;
        for(i_rec = from_rec; pos < n_nucleotides; i_rec++) {
                /* initialize key */
                key.data = &i_rec;
                key.size = sizeof(db_recno_t);

                /* retreive the next rec */
                if ((ret = dbp->get(dbp, NULL, &key, &data, 0)) != 0) {
                        buf[pos] = '\0';
                        return ret;
                }

                while (pos < n_nucleotides && i_rec_offset < data.size)
                {
                        if (is_nucleotide(((char *)data.data)[i_rec_offset])) {
                                buf[pos++] = ((char *)data.data)[i_rec_offset++];
                        }
                        else {
                                i_rec_offset++;
                        }
                }
                /* reset pointer */
                i_rec_offset = 0;
        }

        return 0;
}

// search for pwm match
////////////////////////////////////////////////////////////////////////////////

typedef struct {
        DB* dbp;
        const char* db_name;
        const matrix_t* pwm;
        double threshold;
        int thread_id;
} pthread_pwm_data;

static void* hdb_search_pwm_thread(void* data_)
{
        pthread_pwm_data* data = (pthread_pwm_data*)data_;
        DB* dbp                = data->dbp;
        const char* db_name    = data->db_name;
        const matrix_t* pwm    = data->pwm;
        const double threshold = data->threshold;
        char buf[HUMANDB_RECORD_LENGTH+1];
        size_t pos, i, j, n = 0;

        // get sequence
        for (pos  = 0; hdb_get_sequence(dbp, pos, HUMANDB_RECORD_LENGTH, buf) == 0;
             pos += HUMANDB_RECORD_LENGTH - pwm->columns + 1)
        {
                // loop through sequence
                for (i = 0; i < HUMANDB_RECORD_LENGTH - pwm->columns + 1; i++) {
                        // test pattern
                        double sum = 0;
                        for (j = 0; j < pwm->columns; j++) {
                                if (buf[i+j] == '\0') {
                                        // reached end of upstream sequence or chromosome
                                        printf("%s: Found %d sequences.\n", db_name, (int)n);
                                        return NULL;
                                }
                                else {
                                        switch (buf[i+j]) {
                                        case 'A':
                                        case 'a': sum += pwm->mat[1][j]; break;
                                        case 'C':
                                        case 'c': sum += pwm->mat[3][j]; break;
                                        case 'G':
                                        case 'g': sum += pwm->mat[0][j]; break;
                                        case 'T':
                                        case 't': sum += pwm->mat[2][j]; break;
                                        default: break;
                                        }
                                        if (j == pwm->columns - 1 && sum > threshold) {
                                                printf("(%12s,%010lu,%3.5f): %.*s\n", db_name, (unsigned long)pos+i, sum, pwm->columns, buf+i);
                                                fflush(stdout); n++;
                                        }
                                }
                        }
                }
        }

        return NULL;
}

int hdb_search_pwm(DB* dbp_list[], const int dbp_list_n, const char* db_names[], const matrix_t* pwm, const double threshold)
{
        pthread_t threads[_num_threads];
        pthread_pwm_data data[dbp_list_n];
        int i, j, rc;

        for (i = 0; i < dbp_list_n; i++) {
                data[i].dbp       = dbp_list[i];
                data[i].db_name   = db_names[i];
                data[i].pwm       = pwm;
                data[i].threshold = threshold;
                data[i].thread_id = i;
        }
        for (i = 0; i < dbp_list_n; i += _num_threads) {
                for (j = 0; j < _num_threads && i+j < dbp_list_n; j++) {

                        rc = pthread_create(&threads[j], NULL, hdb_search_pwm_thread, (void *)&data[i+j]);
                        if (rc) {
                                std_err(NONE, "Couldn't create thread.");
                        }
                }
                for (j = 0; j < _num_threads && i+j < dbp_list_n; j++) {
                        rc = pthread_join(threads[j], NULL);
                        if (rc) {
                                std_err(NONE, "Couldn't join thread.");
                        }
                }
        }

        return 0;
}

// search for sequence
////////////////////////////////////////////////////////////////////////////////

typedef struct {
        DB* dbp;
        const char* db_name;
        const char* pattern;
        int pattern_n;
        int thread_id;
} pthread_data;

static inline
char nucleotide_matches(char a, char b) {
        if (tolower(a) == tolower(b)) {
                return 1;
        }
        else {
                return 0;
        }
}

static void* hdb_search_thread(void* data_)
{
        pthread_data* data  = (pthread_data*)data_;
        DB* dbp             = data->dbp;
        const char* db_name = data->db_name;
        const char* pattern = data->pattern;
        const int pattern_n = data->pattern_n;
        char buf[HUMANDB_RECORD_LENGTH+1];
        size_t pos, i, j;

        // get sequence
        for (pos  = 0; hdb_get_sequence(dbp, pos, HUMANDB_RECORD_LENGTH, buf) == 0;
             pos += HUMANDB_RECORD_LENGTH - pattern_n + 1)
        {
                // loop through sequence
                for (i = 0; i < HUMANDB_RECORD_LENGTH - pattern_n + 1; i++) {
                        // test pattern
                        for (j = 0; j < pattern_n; j++) {
                                if (buf[i+j] == '\0') {
                                        // reached end of upstream sequence or chromosome
                                        return NULL;
                                }
                                else if (!nucleotide_matches(buf[i+j], pattern[j])) {
                                        break;
                                }
                                else if (j == pattern_n - 1) {
                                        printf("%s: %010lu\n", db_name, (unsigned long)pos+i);
                                        fflush(stdout);
                                }
                        }
                }
        }

        return NULL;
}

int hdb_search(DB* dbp_list[], const int dbp_list_n, const char* db_names[], const char* pattern, const int pattern_n)
{
        pthread_t threads[_num_threads];
        pthread_data data[dbp_list_n];
        int i, j, rc;

        for (i = 0; i < dbp_list_n; i++) {
                data[i].dbp       = dbp_list[i];
                data[i].db_name   = db_names[i];
                data[i].pattern   = pattern;
                data[i].pattern_n = pattern_n;
                data[i].thread_id = i;
        }

        for (i = 0; i < dbp_list_n; i += _num_threads) {
                for (j = 0; j < _num_threads && i+j < dbp_list_n; j++) {
                        rc = pthread_create(&threads[j], NULL, hdb_search_thread, (void *)&data[i+j]);
                        if (rc) {
                                std_err(NONE, "Couldn't create thread.");
                        }

                }
                for (j = 0; j < _num_threads && i+j < dbp_list_n; j++) {
                        rc = pthread_join(threads[j], NULL);
                        if (rc) {
                                std_err(NONE, "Couldn't join thread.");
                        }
                }
        }

        return 0;
}
