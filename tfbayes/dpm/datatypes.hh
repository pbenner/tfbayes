/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef __TFBAYES_DPM_DATATYPES_HH__
#define __TFBAYES_DPM_DATATYPES_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <vector>
#include <string>

#include <tfbayes/utility/clonable.hh>
#include <tfbayes/utility/linalg.hh>
#include <tfbayes/dpm/index.hh>

// cluster structures
////////////////////////////////////////////////////////////////////////////////

// data types for identifying clusters and baseline models (internal
// use only!)
typedef ssize_t cluster_tag_t;
typedef ssize_t baseline_tag_t;

// a model can be identified by its name and length (this holds for
// baseline and background models)
struct model_id_t {
        std::string name;
        size_t length;

        friend bool operator==(const model_id_t& lhs, const model_id_t& rhs) {
                return lhs.name == rhs.name && lhs.length == rhs.length;
        }
        friend bool operator!=(const model_id_t& lhs, const model_id_t& rhs) {
                return !(lhs == rhs);
        }
};

typedef enum {
        cluster_event_empty, cluster_event_nonempty,
        cluster_event_add_word, cluster_event_remove_word
} cluster_event_t;

#endif /* __TFBAYES_DPM_DATATYPES_HH__ */
