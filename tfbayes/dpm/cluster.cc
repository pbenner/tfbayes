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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <list>

#include <tfbayes/dpm/cluster.hh>

using namespace std;

cluster_t::cluster_t(component_model_t* model, cluster_tag_t cluster_tag, baseline_tag_t baseline_tag,
                 bool destructible, bool record)
        : _model(model), _cluster_tag(cluster_tag), _baseline_tag(baseline_tag),
          _destructible(destructible), _record(record), _size(0)
{ }

cluster_t::cluster_t(component_model_t* model, cluster_tag_t cluster_tag, baseline_tag_t baseline_tag,
                 Observer<cluster_event_t>* observer, bool destructible, bool record)
        : _model(model), _cluster_tag(cluster_tag), _baseline_tag(baseline_tag),
          _destructible(destructible), _record(record), _size(0)
{
        set_observer(observer);
}

cluster_t::cluster_t(const cluster_t& cluster)
        : _model(NULL),
          _cluster_tag(cluster._cluster_tag),
          _baseline_tag(cluster._baseline_tag),
          _destructible(cluster._destructible),
          _record(cluster._record)
{
        operator=(cluster);
}

cluster_t::~cluster_t() {
        delete(_model);
}

cluster_t&
cluster_t::operator=(const cluster_t& cluster)
{
        if (_model) {
                delete(_model);
        }
        _model    = cluster._model->clone();
        _elements = cluster._elements;
        _size     = cluster._size;

        set_observer(cluster.observer);

        return *this;
}

void
cluster_t::add_observations(const range_t& range)
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
cluster_t::remove_observations(const range_t& range)
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

cluster_t::elements_t
cluster_t::elements() const
{
        return _elements;
}

size_t
cluster_t::size() const
{
        return _size;
}

cluster_tag_t
cluster_t::cluster_tag() const
{
        return _cluster_tag;
}

baseline_tag_t
cluster_t::baseline_tag() const
{
        return _baseline_tag;
}

bool
cluster_t::destructible() const
{
        return _destructible;
}

component_model_t&
cluster_t::model()
{
        return *_model;
}

const component_model_t&
cluster_t::model() const
{
        return *_model;
}

ostream&
operator<< (ostream& o, const cluster_t& cluster)
{
        return o << "("
                 << cluster._cluster_tag
                 << ":"
                 << cluster.model().id().name
                 << ":"
                 << cluster.model().id().length
                 << "):"
                 << cluster._size;
}
