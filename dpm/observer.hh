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

#ifndef OBSERVER_HH
#define OBSERVER_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <datatypes.hh>

template <class E> class Observer;
template <class E> class Observed;

template <class E>
class Observer {
public:
        virtual ~Observer() {}
        virtual void update(Observed<E>* observed, E event) = 0;
        virtual void update(Observed<E>* observed, E event, const word_t& word) = 0;
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
        void notify(E event, const word_t& word) {
                if (observer) {
                        observer->update(this, event, word);
                }
        }

        Observer<E>* observer;
};

#endif /* OBSERVER_HH */
