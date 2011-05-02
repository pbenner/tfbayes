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

#include "cluster.hh"

Cluster::Cluster(Data& data)
        : clusters(data.size()),
          assignments(data.size())
{
        for (Cluster::size_type i = 0; i < data.size(); i++) {
                clusters[i].tag = i;
                free_clusters.push_back(&clusters[i]);
        }

        // assign the data to some random clusters
        for (Cluster::size_type i = 0; i < data.size(); i++) {
                Cluster::cluster_tag_t c = rand() % _INIT_NUM_CLASSES;
                Data::element& e = data[i];
                assign(e, c);
        }
}

Cluster::~Cluster() {

}

Cluster::cluster& Cluster::operator[](int c) {
        Cluster::iterator it = begin();
        for (int i = 0; i < c; i++) { it++; }

        return **it;
}

const Cluster::cluster& Cluster::operator[](int c) const {
        Cluster::const_iterator it = begin();
        for (int i = 0; i < c; i++) { it++; }

        return **it;
}

Cluster::cluster_tag_t Cluster::getClusterTag(const Data::element& element) const {
        return assignments[element.tag];
}

void Cluster::release(Data::element& element) {
        Cluster::size_type old_cluster_tag = assignments[element.tag];
        Cluster::elements_t::iterator it = clusters[old_cluster_tag].elements.begin();
        clusters[old_cluster_tag].elements.remove(&element);
        assignments[element.tag] = -1;
        if (clusters[old_cluster_tag].elements.size() == 0) {
                used_clusters.remove(&clusters[old_cluster_tag]);
                free_clusters.push_back(&clusters[old_cluster_tag]);
        }
}

void Cluster::assign(Data::element& element, Cluster::cluster_tag_t cluster_tag) {
        if (clusters[cluster_tag].elements.size() == 0) {
                used_clusters.push_back(&clusters[cluster_tag]);
                free_clusters.remove(&clusters[cluster_tag]);
        }
        assignments[element.tag] = cluster_tag;
        clusters[cluster_tag].elements.push_back(&element);
}

void Cluster::reassign(Data::element& element, Cluster::cluster_tag_t cluster_tag) {
        release(element);
        assign (element, cluster_tag);
}

ostream& operator<< (ostream& o, Cluster const& clusters)
{
        for (Cluster::const_iterator it1 = clusters.begin();
             it1 != clusters.end(); it1++) {
                o << "Cluster: " << (*it1)->tag << endl;
                for (Cluster::elements_t::const_iterator it2 = (*it1)->elements.begin();
                     it2 != (*it1)->elements.end(); it2++) {
                        if (it2 != (*it1)->elements.begin()) {
                                o << ", ";
                        }
                        o << **it2;
                }
                o << endl << endl;
        }

        o << "Used clusters: ";
        for (list<Cluster::cluster*>::const_iterator it = clusters.used_clusters.begin();
             it != clusters.used_clusters.end(); it++) {
                if (it != clusters.used_clusters.begin()) {
                        o << ", ";
                }
                o << (*it)->tag;
        }
        o << endl;

        o << "Free clusters: ";
        for (list<Cluster::cluster*>::const_iterator it = clusters.free_clusters.begin();
             it != clusters.free_clusters.end(); it++) {
                if (it != clusters.free_clusters.begin()) {
                        o << ", ";
                }
                o << (*it)->tag;
        }
        o << endl;

        return o;
}
