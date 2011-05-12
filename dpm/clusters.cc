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

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <list>
using namespace std;

#include "clusters.hh"

Clusters::Clusters(Data& data)
        : clusters(data.size()),
          assignments(data.size()),
          total_elements(0)
{
        // initialize clusters
        for (Clusters::size_type i = 0; i < data.size(); i++) {
                clusters[i].tag = i;
                free_clusters.push_back(&clusters[i]);
        }

        // assign the data to some random clusters
        for (Clusters::size_type i = 0; i < data.size(); i++) {
                Clusters::cluster_tag_t c = rand() % _INIT_NUM_CLASSES;
                Data::element& e = data[i];
                assign(e, c);
        }
}

Clusters::~Clusters() {
}

Clusters::cluster& Clusters::operator[](int c) {
        return clusters[c];
}

const Clusters::cluster& 
Clusters::operator[](int c) const {
        return clusters[c];
}

Clusters::cluster_tag_t
Clusters::getClusterTag(const Data::element& element) const {
        return assignments[element.tag];
}

void
Clusters::release(Data::element& element) {
        Clusters::cluster_tag_t old_cluster_tag = assignments[element.tag];
        if (old_cluster_tag != -1) {
                Clusters::elements_t::iterator it = clusters[old_cluster_tag].elements.begin();
                clusters[old_cluster_tag].elements.remove(&element);
                assignments[element.tag] = -1;
                if (clusters[old_cluster_tag].elements.size() == 0) {
                        used_clusters.remove(&clusters[old_cluster_tag]);
                        free_clusters.push_back(&clusters[old_cluster_tag]);
                }
                total_elements--;
        }
}

void
Clusters::assign(Data::element& element, Clusters::cluster_tag_t cluster_tag) {
        if (clusters[cluster_tag].elements.size() == 0) {
                used_clusters.push_back(&clusters[cluster_tag]);
                free_clusters.remove(&clusters[cluster_tag]);
        }
        assignments[element.tag] = cluster_tag;
        clusters[cluster_tag].elements.push_back(&element);
        total_elements++;
}

void
Clusters::reassign(Data::element& element, Clusters::cluster_tag_t cluster_tag) {
        release(element);
        assign (element, cluster_tag);
}

ostream&
operator<< (ostream& o, Clusters const& clusters)
{
        for (Clusters::const_iterator it1 = clusters.begin();
             it1 != clusters.end(); it1++) {
                o << "Clusters: " << (*it1)->tag << endl;
                for (Clusters::elements_t::const_iterator it2 = (*it1)->elements.begin();
                     it2 != (*it1)->elements.end(); it2++) {
                        if (it2 != (*it1)->elements.begin()) {
                                o << ", ";
                        }
                        o << **it2;
                }
                o << endl << endl;
        }

        o << "Used clusters: ";
        for (list<Clusters::cluster*>::const_iterator it = clusters.used_clusters.begin();
             it != clusters.used_clusters.end(); it++) {
                if (it != clusters.used_clusters.begin()) {
                        o << ", ";
                }
                o << (*it)->tag;
        }
        o << endl;

        o << "Free clusters: ";
        for (list<Clusters::cluster*>::const_iterator it = clusters.free_clusters.begin();
             it != clusters.free_clusters.end(); it++) {
                if (it != clusters.free_clusters.begin()) {
                        o << ", ";
                }
                o << (*it)->tag;
        }
        o << endl;

        return o;
}
