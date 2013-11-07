/* Copyright (C) 2013 Philipp Benner
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

#include <variational.hh>

using namespace std;

// normal factor node
////////////////////////////////////////////////////////////////////////////////

normal_fnode_t::normal_fnode_t(double mean, double precision) :
        dmean(mean),
        dprecision(precision) {
        // set the inbox to the given parameters, however,
        // once a node is connected to a slot, the respective
        // parameter (represented by the dirac distribution)
        // is replaced
        cout << "initializing mailbox" << endl;
        _inbox[1].replace(dmean);
        _inbox[2].replace(dprecision);
        cout << "done." << endl;
}

normal_fnode_t*
normal_fnode_t::clone() const {
        return new normal_fnode_t(*this);
}

bool
normal_fnode_t::link(const std::string& id, variable_node_i& variable_node) {
        if      (id == "output"   ) return base_t::link(0, variable_node);
        else if (id == "mean"     ) return base_t::link(1, variable_node);
        else if (id == "precision") return base_t::link(2, variable_node);
        else return false;
}

bool
normal_fnode_t::is_conjugate(size_t i, variable_node_i& variable_node) const {
        switch (i) {
        case 0: return variable_node.type() == typeid(normal_distribution_t);
        case 1: return variable_node.type() == typeid(normal_distribution_t);
        case 2: return variable_node.type() == typeid( gamma_distribution_t);
        default: assert(false);
        }
}

const p_message_t&
normal_fnode_t::initial_message(size_t i) const {
        switch (i) {
        case 0: return distribution1;
        case 1: return distribution2;
        case 2: return distribution3;
        default: assert(false);
        }
}

const p_message_t&
normal_fnode_t::message(size_t i) {
        switch (i) {
        case 0: return message1();
        case 1: return message2();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
normal_fnode_t::message1() {
        double mean      = _inbox[1]().moment<1>();
        double precision = _inbox[2]().moment<1>();

        cout << "preparing message 1:" << endl
             << "-> mean     : " << mean      << endl
             << "-> precisoin: " << precision << endl;

        distribution1 = normal_distribution_t(mean, precision);

        return distribution1;
}

const p_message_t&
normal_fnode_t::message2() {
        double mean      = _inbox[0]().moment<1>();
        double precision = _inbox[2]().moment<1>();

        cout << "preparing message 2:" << endl
             << "-> mean     : " << mean      << endl
             << "-> precisoin: " << precision << endl;

        distribution2 = normal_distribution_t(mean, precision);

        return distribution2;
}

const p_message_t&
normal_fnode_t::message3() {
        // moments
        double y   = _inbox[0]().moment<1>();
        double y2  = _inbox[0]().moment<2>();
        double mu  = _inbox[1]().moment<1>();
        double mu2 = _inbox[1]().moment<2>();
        // parameters of the gamma distribution
        double shape = 1.5;
        double rate  = 0.5*(y2 - 2.0*y*mu + mu2);

        cout << "preparing message 3:" << endl
             << "-> shape: " << shape << endl
             << "-> rate : " << rate  << endl;

        // replace distribution
        distribution3 = gamma_distribution_t(shape, rate);

        return distribution3;
}

// gamma factor node
////////////////////////////////////////////////////////////////////////////////

gamma_fnode_t::gamma_fnode_t(double shape, double rate) :
        dshape (shape),
        drate  (rate) {
        // set the inbox to the given parameters, however,
        // once a node is connected to a slot, the respective
        // parameter (represented by the dirac distribution)
        // is replaced
        _inbox[1].replace(dshape);
        _inbox[2].replace(drate);
}

gamma_fnode_t*
gamma_fnode_t::clone() const {
        return new gamma_fnode_t(*this);
}

bool
gamma_fnode_t::link(const std::string& id, variable_node_i& variable_node) {
        if      (id == "output") return base_t::link(0, variable_node);
        else if (id == "shape" ) return base_t::link(1, variable_node);
        else if (id == "rate"  ) return base_t::link(2, variable_node);
        else return false;
}

bool
gamma_fnode_t::is_conjugate(size_t i, variable_node_i& variable_node) const {
        switch (i) {
        case 0: return variable_node.type() == typeid(gamma_distribution_t);
        case 1: return false;
        case 2: return variable_node.type() == typeid(gamma_distribution_t);
        default: assert(false);
        }
}

const p_message_t&
gamma_fnode_t::initial_message(size_t i) const {
        switch (i) {
        case 0: return distribution1;
        case 2: return distribution3;
        default: assert(false);
        }
}

const p_message_t&
gamma_fnode_t::message(size_t i) {
        switch (i) {
        case 0: return message1();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
gamma_fnode_t::message1() {
        double shape = _inbox[1]().moment<1>();
        double rate  = _inbox[2]().moment<1>();

        distribution1 = gamma_distribution_t(shape, rate);

        return distribution1;
}

const p_message_t&
gamma_fnode_t::message3() {
        double shape =   _inbox[1]().moment<1>() + 1.0;
        double rate  = - _inbox[0]().moment<1>();

        distribution3 = gamma_distribution_t(shape, rate);

        return distribution3;
}
