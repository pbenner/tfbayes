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
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "data.hh"

Data::Data() {
//        shuffle();
}

Data::~Data() {

}

void Data::shuffle() {
        random_shuffle(elements.begin(), elements.end());
}

const Data::element& Data::operator[](size_t i) const {
        return elements[i];
}

Data::element& Data::operator[](size_t i) {
        return elements[i];
}

ostream& operator<< (ostream& o, Data::element const& element) {
        if (element.x == '\0') {
                o << "\n";
        }
        else {
                o << element.tag << ":" << element.x;
                o << "@" << element.original_cluster;
        }

        return o;
}

ostream& operator<< (ostream& o, Data const& data) {
        for (Data::const_iterator it = data.begin(); it != data.end(); it++) {
                if (it != data.begin()) {
                        o << ", ";
                }
                o << *it;
        }

        return o;
}
