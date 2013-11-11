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

#include <tfbayes/fg/variational.hh>
#include <tfbayes/utility/debug.hh>

using namespace std;

// normal factor node
////////////////////////////////////////////////////////////////////////////////

normal_fnode_t::normal_fnode_t(double mean, double precision,
                               const std::string& name) :
        base_t       (name),
        dmean        (mean),
        dprecision   (precision),
        // initial distributions
        distribution1(0,0.01),
        distribution2(0,0.01),
        distribution3(1,0.5) {
        // set the inbox to the given parameters, however,
        // once a node is connected to a slot, the respective
        // parameter (represented by the dirac distribution)
        // is replaced
        _inbox[1].replace(dmean);
        _inbox[2].replace(dprecision);
}

normal_fnode_t::normal_fnode_t(const normal_fnode_t& normal_fnode) :
        base_t       (normal_fnode),
        dmean        (normal_fnode.dmean),
        dprecision   (normal_fnode.dprecision),
        distribution1(normal_fnode.distribution1),
        distribution2(normal_fnode.distribution2),
        distribution3(normal_fnode.distribution3) {
        _inbox[1].replace(dmean);
        _inbox[2].replace(dprecision);
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
        assert(_inbox[0]().dimension() == 1);
        assert(_inbox[1]().dimension() == 1);
        assert(_inbox[2]().dimension() == 1);
        switch (i) {
        case 0: return message1();
        case 1: return message2();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
normal_fnode_t::message1() {
        double mean      = _inbox[1]().moment<1>()[0];
        double precision = _inbox[2]().moment<1>()[0];

        debug("normal message 1 (normal): " << this->name() << endl);
        distribution1 = normal_distribution_t(mean, precision);
        debug(endl);

        return distribution1;
}

const p_message_t&
normal_fnode_t::message2() {
        double mean      = _inbox[0]().moment<1>()[0];
        double precision = _inbox[2]().moment<1>()[0];

        debug("normal message 2 (normal): " << this->name() << endl);
        distribution2 = normal_distribution_t(mean, precision);
        debug(endl);

        return distribution2;
}

const p_message_t&
normal_fnode_t::message3() {
        // moments
        double y   = _inbox[0]().moment<1>()[0];
        double y2  = _inbox[0]().moment<2>()[0];
        double mu  = _inbox[1]().moment<1>()[0];
        double mu2 = _inbox[1]().moment<2>()[0];
        // parameters of the gamma distribution
        double shape = 1.5;
        double rate  = 0.5*(y2 - 2.0*y*mu + mu2);
        assert(rate > 0.0);

        debug("normal message 3 (gamma): " << this->name() << endl);
        debug("-> y    : " << y     << endl);
        debug("-> y2   : " << y2    << endl);
        debug("-> mu   : " << mu    << endl);
        debug("-> mu2  : " << mu2   << endl);
        debug("-> shape: " << shape << endl);
        debug("-> rate : " << rate  << endl);
        // replace distribution
        distribution3 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution3;
}

// normal factor node
////////////////////////////////////////////////////////////////////////////////

pnormal_fnode_t::pnormal_fnode_t(size_t dim, double mean, double precision,
                                 const std::string& name) :
        base_t       (name),
        dmean        (mean),
        dprecision   (precision),
        // initial distributions
        distribution1(dim,0,0.01),
        distribution2(0,0.01),
        distribution3(1,0.5),
        dimension    (dim) {
        assert(dim > 0);
        // set the inbox to the given parameters, however,
        // once a node is connected to a slot, the respective
        // parameter (represented by the dirac distribution)
        // is replaced
        _inbox[1].replace(dmean);
        _inbox[2].replace(dprecision);
}

pnormal_fnode_t::pnormal_fnode_t(const pnormal_fnode_t& pnormal_fnode) :
        base_t       (pnormal_fnode),
        dmean        (pnormal_fnode.dmean),
        dprecision   (pnormal_fnode.dprecision),
        distribution1(pnormal_fnode.distribution1),
        distribution2(pnormal_fnode.distribution2),
        distribution3(pnormal_fnode.distribution3),
        dimension    (pnormal_fnode.dimension) {
        _inbox[1].replace(dmean);
        _inbox[2].replace(dprecision);
}

pnormal_fnode_t*
pnormal_fnode_t::clone() const {
        return new pnormal_fnode_t(*this);
}

bool
pnormal_fnode_t::link(const std::string& id, variable_node_i& variable_node) {
        if      (id == "output"   ) return base_t::link(0, variable_node);
        else if (id == "mean"     ) return base_t::link(1, variable_node);
        else if (id == "precision") return base_t::link(2, variable_node);
        else return false;
}

bool
pnormal_fnode_t::is_conjugate(size_t i, variable_node_i& variable_node) const {
        switch (i) {
        case 0: return variable_node.type() == typeid(pnormal_distribution_t);
        case 1: return variable_node.type() == typeid( normal_distribution_t);
        case 2: return variable_node.type() == typeid(  gamma_distribution_t);
        default: assert(false);
        }
}

const p_message_t&
pnormal_fnode_t::initial_message(size_t i) const {
        switch (i) {
        case 0: return distribution1;
        case 1: return distribution2;
        case 2: return distribution3;
        default: assert(false);
        }
}

const p_message_t&
pnormal_fnode_t::message(size_t i) {
        assert(_inbox[0]().dimension() == dimension);
        assert(_inbox[1]().dimension() == 1);
        assert(_inbox[2]().dimension() == 1);
        switch (i) {
        case 0: return message1();
        case 1: return message2();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
pnormal_fnode_t::message1() {
        double mean      = _inbox[1]().moment<1>()[0];
        double precision = _inbox[2]().moment<1>()[0];

        debug("normal message 1 (normal): " << this->name() << endl);
        distribution1 = pnormal_distribution_t(dimension, mean, precision);
        debug(endl);

        return distribution1;
}

const p_message_t&
pnormal_fnode_t::message2() {
        double d         = static_cast<double>(dimension);
        double mean      = 0.0;
        for (size_t i = 0; i < dimension; i++) {
                mean += 1.0/d*_inbox[0]().moment<1>()[i];
        }
        double precision = d*_inbox[2]().moment<1>()[0];

        debug("normal message 2 (normal): " << this->name() << endl);
        distribution2 = normal_distribution_t(mean, precision);
        debug(endl);

        return distribution2;
}

const p_message_t&
pnormal_fnode_t::message3() {
        // moments
        double d   = static_cast<double>(dimension);
        double y   = 0.0;
        double y2  = 0.0;
        for (size_t i = 0; i < dimension; i++) {
                y  += _inbox[0]().moment<1>()[i];
                y2 += _inbox[0]().moment<2>()[i];
        }
        double mu  = _inbox[1]().moment<1>()[0];
        double mu2 = _inbox[1]().moment<2>()[0];
        // parameters of the gamma distribution
        double shape = d/2.0 + 1.0;
        double rate  = 0.5*(y2 - 2.0*y*mu + d*mu2);
        assert(rate > 0.0);

        debug("normal message 3 (gamma): " << this->name() << endl);
        debug("-> dim  : " << d     << endl);
        debug("-> y    : " << y     << endl);
        debug("-> y2   : " << y2    << endl);
        debug("-> mu   : " << mu    << endl);
        debug("-> mu2  : " << mu2   << endl);
        debug("-> shape: " << shape << endl);
        debug("-> rate : " << rate  << endl);
        // replace distribution
        distribution3 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution3;
}

// gamma factor node
////////////////////////////////////////////////////////////////////////////////

gamma_fnode_t::gamma_fnode_t(double shape, double rate,
                             const std::string& name) :
        base_t (name),
        dshape (shape),
        drate  (rate),
        // initial distributions
        distribution1(1,1),
        distribution3(1,1) {
        // set the inbox to the given parameters, however,
        // once a node is connected to a slot, the respective
        // parameter (represented by the dirac distribution)
        // is replaced
        _inbox[1].replace(dshape);
        _inbox[2].replace(drate);
}

gamma_fnode_t::gamma_fnode_t(const gamma_fnode_t& gamma_fnode) :
        base_t       (gamma_fnode),
        dshape       (gamma_fnode.dshape),
        drate        (gamma_fnode.drate),
        distribution1(gamma_fnode.distribution1),
        distribution3(gamma_fnode.distribution3) {
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
        assert(_inbox[0]().dimension() == 1);
        assert(_inbox[1]().dimension() == 1);
        assert(_inbox[2]().dimension() == 1);
        switch (i) {
        case 0: return message1();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
gamma_fnode_t::message1() {
        double shape = _inbox[1]().moment<1>()[0];
        double rate  = _inbox[2]().moment<1>()[0];

        debug("gamma message 1 (gamma): " << this->name() << endl);
        distribution1 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution1;
}

const p_message_t&
gamma_fnode_t::message3() {
        double shape =   _inbox[1]().moment<1>()[0] + 1.0;
        double rate  = - _inbox[0]().moment<1>()[0];

        debug("gamma message 1 (gamma): " << this->name() << endl);
        distribution3 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution3;
}
