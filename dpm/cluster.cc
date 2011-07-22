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

Cluster::Cluster(ComponentModel* model, cluster_tag_t tag)
        : _model(model), _tag(tag), _destructible(true), _size(0)
{ }

Cluster::Cluster(ComponentModel* model, cluster_tag_t tag, bool destructible)
        : _model(model), _tag(tag), _destructible(destructible), _size(0)
{ }

Cluster::Cluster(ComponentModel* model, cluster_tag_t tag, Observer<cluster_event_t>* observer)
        : _model(model), _tag(tag), _destructible(true), _size(0)
{
        set_observer(observer);
}

Cluster::Cluster(ComponentModel* model, cluster_tag_t tag, Observer<cluster_event_t>* observer, bool destructible)
        : _model(model), _tag(tag), _destructible(destructible), _size(0)
{
        set_observer(observer);
}

Cluster::Cluster(const Cluster& cluster)
        : _model(cluster._model->clone()),
          _tag(cluster._tag),
          _destructible(cluster._destructible),
          _size(cluster._size)
{
        set_observer(cluster.observer);
}

Cluster::~Cluster() {
        delete(_model);
}

void
Cluster::add_observations(const range_t& range)
{
        if (_size == 0) {
                _size += _model->add(range);
                notify(cluster_event_nonempty);
        }
        else {
                _size += _model->add(range);
        }
        notify(cluster_event_add_word, range);
}

void
Cluster::remove_observations(const range_t& range)
{
        if (_size >= _model->count(range)) {
                _size -= _model->remove(range);
                if (_size == 0) {
                        notify(cluster_event_empty);
                }
                notify(cluster_event_remove_word, range);
        }
}

size_t
Cluster::size() const
{
        return _size;
}

cluster_tag_t
Cluster::tag() const
{
        return _tag;
}

bool
Cluster::destructible() const
{
        return _destructible;
}

const ComponentModel&
Cluster::model() const
{
        return *_model;
}

ostream&
operator<< (ostream& o, const Cluster& cluster)
{
        return o << cluster._tag << "@" << cluster._size;
}
