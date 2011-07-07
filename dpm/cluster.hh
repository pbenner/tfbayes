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
#include "statistics.hh"

////////////////////////////////////////////////////////////////////////////////

typedef size_t cluster_tag_t;

typedef enum {
        cluster_event_empty, cluster_event_nonempty
} cluster_event_t;

template <class E> class Observer;
template <class E> class Observed;

template <class E>
class Observer {
public:
        virtual void update(Observed<E>* observed, E event) = 0;
};

template <class E>
class Observed {
public:
        Observed() : observer(NULL) {};

        void set_observer(Observer<E>* observer) {
                this->observer = observer;
        }

protected:
        void notify(E event) {
                if (observer) {
                        observer->update(this, event);
                }
        }

private:
        Observer<E>* observer;
};

// This class represents a single cluster for the dirichlet process
// mixture. Each cluster is linked to a probability distribution and
// keeps track of the sufficient statistics assiciated to each cluster.
////////////////////////////////////////////////////////////////////////////////

class Cluster : public Observed<cluster_event_t> {
public:
         Cluster(Distribution* distribution, cluster_tag_t tag);
         Cluster(Distribution* distribution, cluster_tag_t tag, bool destructible);
         Cluster(Distribution* distribution, cluster_tag_t tag, Observer<cluster_event_t>* observer);
         Cluster(Distribution* distribution, cluster_tag_t tag, bool destructible, Observer<cluster_event_t>* observer);
        ~Cluster();

        // methods
        void add_word(const word_t& word);
        void remove_word(const word_t& word);

        size_t size() const;

        // public variables
        Distribution* distribution;
        cluster_tag_t tag;
        // is this cluster essential for the mixture?
        // i.e. for the background model
        bool destructible;
private:
        size_t _size;
};

#endif /* CLUSTER_HH */
