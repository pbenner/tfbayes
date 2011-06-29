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
#include <math.h>

#include <pthread.h>
#include <db.h>

#include <tfbayes/exception.h>

#include "humandb.h"

static
void add_codon(char codon[], long result[])
{
        int i, pos = 0;

        for (i = 0; i < 3; i++) {
                switch (codon[i]) {
                case 'a':
                case 'A':
                        pos += 0 * pow(4, i); break;
                case 'c':
                case 'C':
                        pos += 1 * pow(4, i); break;
                case 'g':
                case 'G':
                        pos += 2 * pow(4, i); break;
                case 't':
                case 'T':
                        pos += 3 * pow(4, i); break;
                default:
                        return;
                }
        }
        result[pos] += 1;
}


void hdb_count_codons(DB *dbp, long positions[], size_t n, long result[])
{
        int i, j;

        for (i = 0; i*2 < n; i+= 2) {
                long from = positions[i];
                long to   = positions[i+1];
                char buf[to-from+1];
                hdb_get_sequence(dbp, from, to-from, buf);

                for (j = 0; buf[j]; j++) {
                        add_codon(buf+j, result);
                }
        }
}

void hdb_count_codons_upstream(DB *dbp, size_t len, long result[])
{
        int j;
        char buf[len+1];

//        printf("Hello 1\n");
        fflush(stdout);
        hdb_get_sequence_pure(dbp, 0, len, buf);
//        printf("Hello 2\n");
        fflush(stdout);

        for (j = 0; buf[j]; j++) {
                add_codon(buf+j, result);
        }
}
