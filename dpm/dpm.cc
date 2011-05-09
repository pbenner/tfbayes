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

#include "cluster.hh"
#include "data.hh"
#include "dpm.hh"
#include "statistics.hh"

using namespace std;

DPM::DPM(Data* data) : da(data), cl(*da) {
        hist_switches.push_back(0);
        hist_likelihood.push_back(likelihood());
}

DPM::~DPM() {
        delete(da);
}

ostream& operator<< (ostream& o, DPM const& dpm)
{
        return o << dpm.cl;
}
