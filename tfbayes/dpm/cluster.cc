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

#include <boost/format.hpp>

using namespace std;

cluster_t::cluster_t(component_model_t* model, cluster_tag_t cluster_tag, baseline_tag_t baseline_tag,
                     bool destructible, bool record)
        : base_t         ()
        , m_model        (model)
        , m_cluster_tag  (cluster_tag)
        , m_baseline_tag (baseline_tag)
        , m_destructible (destructible)
        , m_record       (record)
        , m_elements     ()
        , m_size         (0)
{ }

cluster_t::cluster_t(component_model_t* model, cluster_tag_t cluster_tag, baseline_tag_t baseline_tag,
                     Observer<cluster_event_t>* observer, bool destructible, bool record)
        : base_t         ()
        , m_model        (model)
        , m_cluster_tag  (cluster_tag)
        , m_baseline_tag (baseline_tag)
        , m_destructible (destructible)
        , m_record       (record)
        , m_elements     ()
        , m_size         (0)
{
        set_observer(observer);
}

cluster_t::cluster_t(const cluster_t& cluster)
        : base_t         (cluster)
        , m_model        (cluster.m_model->clone())
        , m_cluster_tag  (cluster.m_cluster_tag)
        , m_baseline_tag (cluster.m_baseline_tag)
        , m_destructible (cluster.m_destructible)
        , m_record       (cluster.m_record)
        , m_elements     (cluster.m_elements)
        , m_size         (cluster.m_size)
{
        set_observer(cluster.observer);
}

cluster_t::~cluster_t() {
        if (m_model) {
                delete(m_model);
        }
}

cluster_t&
cluster_t::operator=(const cluster_t& cluster)
{
        if (m_model) {
                delete(m_model);
        }
        m_model        = cluster.m_model->clone();
        m_cluster_tag  = cluster.m_cluster_tag;
        m_baseline_tag = cluster.m_baseline_tag;
        m_destructible = cluster.m_destructible;
        m_record       = cluster.m_record;
        m_elements     = cluster.m_elements;
        m_size         = cluster.m_size;

        set_observer(cluster.observer);

        return *this;
}

void
cluster_t::add_observations(const range_t& range)
{
        if (m_size == 0) {
                m_size += m_model->add(range);
                notify(cluster_event_nonempty);
        }
        else {
                m_size += m_model->add(range);
        }
        notify(cluster_event_add_word, range);

        if (m_record) {
                m_elements.insert(range);
        }
}

void
cluster_t::remove_observations(const range_t& range)
{
        if (m_size >= m_model->count(range)) {
                m_size -= m_model->remove(range);
                if (m_size == 0) {
                        notify(cluster_event_empty);
                }
                notify(cluster_event_remove_word, range);

                if (m_record) {
                        m_elements.erase(range);
                }
        }
}

cluster_t::elements_t
cluster_t::elements() const
{
        return m_elements;
}

size_t
cluster_t::size() const
{
        return m_size;
}

cluster_tag_t
cluster_t::cluster_tag() const
{
        return m_cluster_tag;
}

baseline_tag_t
cluster_t::baseline_tag() const
{
        return m_baseline_tag;
}

bool
cluster_t::destructible() const
{
        return m_destructible;
}

component_model_t&
cluster_t::model()
{
        return *m_model;
}

const component_model_t&
cluster_t::model() const
{
        return *m_model;
}

ostream&
operator<< (ostream& o, const cluster_t& cluster)
{
        o << boost::format("Cluster %d:\n"         ) % cluster.cluster_tag()
          << boost::format(" -> model name  : %s\n") % cluster.model().id().name
          << boost::format(" -> model length: %d\n") % cluster.model().id().length
          << boost::format(" -> elements    : ");
        for (cluster_t::iterator it = cluster.begin(); it != cluster.end(); it++) {
                o << *it << " ";
        }
        return o << endl;
}
