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

#ifndef DATATYPES_HH
#define DATATYPES_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

// efficient representation of words (nucleotide sequences)
// of variable length
typedef struct {
        size_t sequence;
        size_t position;
        size_t length;
        const std::vector<std::vector<char> >& sequences;
} word_t;

typedef struct {
        size_t sequence;
        size_t position;
} element_t;

typedef size_t cluster_tag_t;

typedef enum {
        cluster_event_empty, cluster_event_nonempty,
        cluster_event_add_word, cluster_event_remove_word
} cluster_event_t;

typedef struct {
        std::vector<double> switches;
        std::vector<double> likelihood;
        std::vector<size_t> components;
} sampling_history_t;

static inline
char code_nucleotide(char a)
{
        switch (a) {
        case 'A':
        case 'a':
                return 0;
        break;
        case 'C':
        case 'c':
                return 1;
        break;
        case 'G':
        case 'g':
                return 2;
        break;
        case 'T':
        case 't':
                return 3;
        break;
        }
        std::cerr << "code_nucleotide(): found non-nucleotide.\n" << std::endl;
        exit(EXIT_FAILURE);
}

static inline
char decode_nucleotide(char a)
{
        switch (a) {
        case 0:
                return 'a';
        break;
        case 1:
                return 'c';
        break;
        case 2:
                return 'g';
        break;
        case 3:
                return 't';
        break;
        }
        std::cerr << "decode_nucleotide(): found non-nucleotide.\n" << std::endl;
        exit(EXIT_FAILURE);
}

static inline
bool is_nucleotide(char S) {
        if (S == 'A' || S == 'C' || S == 'G' || S == 'T' ||
            S == 'a' || S == 'c' || S == 'g' || S == 't') {
                return true;
        }
        else {
                return false;
        }
}

#endif /* DATATYPES_HH */
