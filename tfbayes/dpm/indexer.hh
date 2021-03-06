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

#ifndef __TFBAYES_DPM_INDEXER_HH__
#define __TFBAYES_DPM_INDEXER_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>

#include <tfbayes/dpm/index.hh>

//
// Interface for (randomly) indexing data structures
//

class indexer_t {
public:
        virtual ~indexer_t() {}

        // type definitions
        typedef std::vector<index_t>::iterator iterator;
        typedef std::vector<index_t>::const_iterator const_iterator;

        typedef std::vector<index_t>::const_iterator sampling_iterator;

        // iterators
        ////////////////////////////////////////////////////////////////////////
        virtual iterator begin() = 0;
        virtual iterator end()   = 0;

        virtual const_iterator begin() const = 0;
        virtual const_iterator end()   const = 0;

        virtual sampling_iterator sampling_begin() const = 0;
        virtual sampling_iterator sampling_end()   const = 0;

        // operators
        ////////////////////////////////////////////////////////////////////////
        virtual size_t elements() const = 0;

        virtual void shuffle() = 0;
};

#endif /* __TFBAYES_DPM_INDEXER_HH__ */
