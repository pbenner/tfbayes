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
        std::vector<std::vector<double> > switches;
        std::vector<std::vector<double> > likelihood;
        std::vector<std::vector<size_t> > components;
} sampling_history_t;

typedef std::vector<std::vector<double> > posterior_t;

#endif /* DATATYPES_HH */
