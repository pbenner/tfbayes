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

#ifndef CLUSTER_HH
#define CLUSTER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <list>
#include <vector>
using namespace std;

#include "data.hh"

class Cluster {
public:
        // constructors and destructors
        Cluster(Data& data);
        ~Cluster();

        // type definitions
        typedef vector<Data::element *>::size_type size_type;
        typedef list<Data::element *> elements_t;
        typedef size_type cluster_tag_t;
        typedef struct {
                elements_t elements;
                cluster_tag_t tag;
        } cluster; 

        typedef list<cluster*>::iterator iterator;
        typedef list<cluster*>::const_iterator const_iterator;

        // iterators
        iterator begin() { return used_clusters.begin(); }
        iterator end()   { return used_clusters.end(); }

        const_iterator begin() const { return used_clusters.begin(); }
        const_iterator end()   const { return used_clusters.end(); }

        // operators
              cluster& operator[](int c);
        const cluster& operator[](int c) const;

        friend ostream& operator<< (std::ostream& o, Cluster const& cluster);

        // methods
        cluster_tag_t getClusterTag(const Data::element& element) const;
        void release (Data::element& element);
        void assign  (Data::element& element, cluster_tag_t cluster_tag);
        void reassign(Data::element& element, cluster_tag_t cluster_tag);

        cluster* next_free_cluster() {
                if (free_clusters.size() >= 1) {
                        return free_clusters.front();
                }
                else {
                        return NULL;
                }
        }

        size_type size() { return used_clusters.size(); }

private:
        static const int _INIT_NUM_CLASSES = 2;

        vector<cluster> clusters;
        vector<size_type> assignments;
        list<cluster*> used_clusters;
        list<cluster*> free_clusters;
};

#endif /* CLUSTER_HH */
