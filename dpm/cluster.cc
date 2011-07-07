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

#include <cluster.hh>

using namespace std;

Cluster::Cluster(Distribution* distribution, cluster_tag_t tag)
        : distribution(distribution), tag(tag), destructible(true), _size(0)
{ }

Cluster::Cluster(Distribution* distribution, cluster_tag_t tag, bool destructible)
        : distribution(distribution), tag(tag), destructible(destructible), _size(0)
{ }

Cluster::Cluster(Distribution* distribution, cluster_tag_t tag, Observer<cluster_event_t>* observer)
        : distribution(distribution), tag(tag), destructible(true), _size(0)
{
        set_observer(observer);
}

Cluster::Cluster(Distribution* distribution, cluster_tag_t tag, bool destructible, Observer<cluster_event_t>* observer)
        : distribution(distribution), tag(tag), destructible(destructible), _size(0)
{
        set_observer(observer);
}

Cluster::~Cluster() {
        delete(distribution);
}

void
Cluster::add_word(const word_t& word)
{
        if (_size == 0) {
                _size += distribution->add_observations(word);
                notify(cluster_event_nonempty);
        }
        else {
                _size += distribution->add_observations(word);
        }
}

void
Cluster::remove_word(const word_t& word)
{
        _size -= distribution->remove_observations(word);
        if (_size == 0) {
                notify(cluster_event_empty);
        }
}

size_t
Cluster::size() const
{
        return _size;
}
