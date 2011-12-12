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

Cluster::Cluster(ComponentModel* model, cluster_tag_t cluster_tag, model_tag_t model_tag,
                 bool destructible, bool record)
        : _model(model), _cluster_tag(cluster_tag), _model_tag(model_tag),
          _destructible(destructible), _record(record), _size(0)
{ }

Cluster::Cluster(ComponentModel* model, cluster_tag_t cluster_tag, model_tag_t model_tag,
                 Observer<cluster_event_t>* observer, bool destructible, bool record)
        : _model(model), _cluster_tag(cluster_tag), _model_tag(model_tag),
          _destructible(destructible), _record(record), _size(0)
{
        set_observer(observer);
}

Cluster::Cluster(const Cluster& cluster)
        : _model(cluster._model->clone()),
          _cluster_tag(cluster._cluster_tag),
          _model_tag(cluster._model_tag),
          _destructible(cluster._destructible),
          _record(cluster._record),
          _size(cluster._size)
{
        set_observer(cluster.observer);
}

Cluster::~Cluster() {
        delete(_model);
}

void
Cluster::operator=(const Cluster& cluster)
{
        delete(_model);
        _model = cluster._model->clone();
        _elements = cluster._elements;
        _size = cluster._size;
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

        if (_record) {
                _elements.insert(range);
        }
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

                if (_record) {
                        _elements.erase(range);
                }
        }
}

Cluster::elements_t
Cluster::elements() const
{
        return _elements;
}

size_t
Cluster::size() const
{
        return _size;
}

cluster_tag_t
Cluster::cluster_tag() const
{
        return _cluster_tag;
}

model_tag_t
Cluster::model_tag() const
{
        return _model_tag;
}

bool
Cluster::destructible() const
{
        return _destructible;
}

ComponentModel&
Cluster::model()
{
        return *_model;
}

ostream&
operator<< (ostream& o, const Cluster& cluster)
{
        return o << "("
                 << cluster._cluster_tag
                 << ":"
                 << cluster._model_tag
                 << "):"
                 << cluster._size;
}
