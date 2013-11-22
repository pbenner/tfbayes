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

#define _USE_MATH_DEFINES
#include <cmath>

#include <tfbayes/fg/variational.hh>
#include <tfbayes/utility/debug.hh>

using namespace std;

// normal factor node
////////////////////////////////////////////////////////////////////////////////

normal_fnode_t::normal_fnode_t(const std::string& name,
                               double mean, double precision, size_t dimension) :
        base_t       (name),
        dmean        (),
        dprecision   (),
        // initial distributions
        distribution1(),
        distribution2(),
        distribution3(),
        dimension    (dimension) {
        // initialize parameters
        dmean      = distribution2.statistics(mean);
        dprecision = distribution3.statistics(precision);
        // initialize links
        _links[1] = boost::bind(&normal_fnode_t::dmean, this);
        _links[2] = boost::bind(&normal_fnode_t::dprecision, this);
        assert(dimension > 0);
        assert(precision > 0.0);
}

normal_fnode_t::normal_fnode_t(const normal_fnode_t& normal_fnode) :
        base_t       (normal_fnode),
        dmean        (normal_fnode.dmean),
        dprecision   (normal_fnode.dprecision),
        distribution1(normal_fnode.distribution1),
        distribution2(normal_fnode.distribution2),
        distribution3(normal_fnode.distribution3),
        dimension    (normal_fnode.dimension) {
        assert(dimension > 0);
        _links[1] = boost::bind(&normal_fnode_t::dmean, this);
        _links[2] = boost::bind(&normal_fnode_t::dprecision, this);
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

double
normal_fnode_t::free_energy() const
{
        double d       = static_cast<double>(dimension);
        double y       = _links[0]()[0];
        double y2      = _links[0]()[1];
        double mu      = _links[1]()[0];
        double mu2     = _links[1]()[1];
        double log_tau = _links[2]()[0];
        double tau     = _links[2]()[1];
        double result  = 0.0;

        // log base measure
        result -= d/2.0*std::log(2.0*M_PI);
        // log partition
        result -= d/2.0*(mu2*tau - log_tau);
        // parameters * statistics
        result +=  mu*tau*y;
        result -= 0.5*tau*y2;

        debug(boost::format("factor node %s:%x computed free energy: %d\n")
              % base_t::name() % this % result);

        return result;
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
normal_fnode_t::operator()(size_t i) {
        switch (i) {
        case 0: return message1();
        case 1: return message2();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
normal_fnode_t::message1() {
        double mean      = _links[1]()[0];
        double precision = _links[2]()[1];

        debug("normal message 1 (normal): " << this->name() << endl);
        distribution1 = normal_distribution_t(mean, precision);
        debug(endl);

        return distribution1;
}

const p_message_t&
normal_fnode_t::message2() {
        double d         = static_cast<double>(dimension);
        double mean      = 1.0/d * _links[0]()[0];
        double precision =     d * _links[2]()[1];

        debug("normal message 2 (normal): " << this->name() << endl);
        distribution2 = normal_distribution_t(mean, precision);
        debug(endl);

        return distribution2;
}

const p_message_t&
normal_fnode_t::message3() {
        // moments
        double d   = static_cast<double>(dimension);
        double y   = _links[0]()[0];
        double y2  = _links[0]()[1];
        double mu  = _links[1]()[0];
        double mu2 = _links[1]()[1];
        // parameters of the gamma distribution
        double shape = d/2.0 + 1.0;
        double rate  = 0.5*(y2 - 2.0*y*mu + d*mu2);
        assert(rate >= 0.0);
        // in case we have some strange data
        rate = std::max(1e-20, rate);

        debug("normal message 3 (gamma): " << this->name() << endl);
        // replace distribution
        distribution3 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution3;
}

// gamma factor node
////////////////////////////////////////////////////////////////////////////////

gamma_fnode_t::gamma_fnode_t(const std::string& name,
                             double shape, double rate, size_t dimension) :
        base_t       (name),
        dshape       (1, shape),
        drate        (),
        // initial distributions
        distribution1(),
        distribution3(),
        dimension    (dimension) {
        assert(dimension > 0);
        // initialize parameters
        drate = distribution3.statistics(rate);
        // initialize links
        _links[1] = boost::bind(&gamma_fnode_t::dshape, this);
        _links[2] = boost::bind(&gamma_fnode_t::drate,  this);
}

gamma_fnode_t::gamma_fnode_t(const gamma_fnode_t& gamma_fnode) :
        base_t       (gamma_fnode),
        dshape       (gamma_fnode.dshape),
        drate        (gamma_fnode.drate),
        distribution1(gamma_fnode.distribution1),
        distribution3(gamma_fnode.distribution3),
        dimension    (gamma_fnode.dimension) {
        assert(dimension > 0);
        _links[1] = boost::bind(&gamma_fnode_t::dshape, this);
        _links[2] = boost::bind(&gamma_fnode_t::drate,  this);
}

gamma_fnode_t*
gamma_fnode_t::clone() const {
        return new gamma_fnode_t(*this);
}

bool
gamma_fnode_t::link(const std::string& id, variable_node_i& variable_node) {
        if      (id == "output") return base_t::link(0, variable_node);
        else if (id == "rate"  ) return base_t::link(2, variable_node);
        else return false;
}

double
gamma_fnode_t::free_energy() const
{
        double d       = static_cast<double>(dimension);
        double log_tau = _links[0]()[0];
        double tau     = _links[0]()[1];
        double shape   = _links[1]()[0];
        double rate    = _links[2]()[1];
        double result  = 0.0;

        // log partition
        result -= d*(boost::math::lgamma(shape) - shape*std::log(rate));
        // parameters * statistics
        result += (shape-1.0)*log_tau;
        result -= rate*tau;

        debug(boost::format("factor node %s:%x computed free energy: %d\n")
              % base_t::name() % this % result);

        return result;
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
gamma_fnode_t::operator()(size_t i) {
        switch (i) {
        case 0: return message1();
        case 2: return message3();
        default: assert(false);
        }
}

const p_message_t&
gamma_fnode_t::message1() {
        double shape = _links[1]()[0];
        double rate  = _links[2]()[1];

        debug("gamma message 1 (gamma): " << this->name() << endl);
        distribution1 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution1;
}

const p_message_t&
gamma_fnode_t::message3() {
        double shape =   _links[1]()[0] + 1.0;
        double rate  = - _links[0]()[1];

        debug("gamma message 1 (gamma): " << this->name() << endl);
        distribution3 = gamma_distribution_t(shape, rate);
        debug(endl);

        return distribution3;
}

// dirichlet factor node
////////////////////////////////////////////////////////////////////////////////

dirichlet_fnode_t::dirichlet_fnode_t(const std::string& name,
                                     const vector_t& alpha,
                                     size_t dimension) :
        base_t       (name),
        dalpha       (alpha),
        // initial distributions
        distribution1(alpha.size()),
        dimension    (dimension) {
        assert(dimension > 0);
}

dirichlet_fnode_t::dirichlet_fnode_t(const dirichlet_fnode_t& dirichlet_fnode) :
        base_t       (dirichlet_fnode),
        dalpha       (dirichlet_fnode.dalpha),
        distribution1(dirichlet_fnode.distribution1),
        dimension    (dirichlet_fnode.dimension) {
        assert(dimension > 0);
}

dirichlet_fnode_t*
dirichlet_fnode_t::clone() const {
        return new dirichlet_fnode_t(*this);
}

bool
dirichlet_fnode_t::link(const std::string& id, variable_node_i& variable_node) {
        if      (id == "output") return base_t::link(0, variable_node);
        else return false;
}

double
dirichlet_fnode_t::free_energy() const
{
        vector_t log_theta = _links[0]();
        assert(log_theta.size() == dalpha.size());
        double result = 0.0;

        // log partition
        double sum = std::accumulate(dalpha.begin(), dalpha.end(), 0.0);
        for (size_t i = 0; i < dalpha.size(); i++) {
                result -= boost::math::lgamma(dalpha[i]);
        }
        result += boost::math::lgamma(sum);
        // statistics & parameters
        for (size_t i = 0; i < dalpha.size(); i++) {
                result += (dalpha[i]-1.0)*log_theta[i];
        }

        debug(boost::format("factor node %s:%x computed free energy: %d\n")
              % base_t::name() % this % result);

        return result;
}

bool
dirichlet_fnode_t::is_conjugate(size_t i, variable_node_i& variable_node) const {
        switch (i) {
        case 0: return variable_node.type() == typeid(dirichlet_distribution_t);
        default: assert(false);
        }
}

const p_message_t&
dirichlet_fnode_t::operator()(size_t i) {
        switch (i) {
        case 0: return message1();
        default: assert(false);
        }
}

const p_message_t&
dirichlet_fnode_t::message1() {
        debug("dirichlet message 1 (dirichlet): " << this->name() << endl);
        debug(endl);

        return distribution1;
}

// categorical factor node
////////////////////////////////////////////////////////////////////////////////

categorical_fnode_t::categorical_fnode_t(const std::string& name,
                                         const vector_t& theta,
                                         size_t dimension) :
        base_t       (name),
        dtheta       (),
        // initial distributions
        distribution1(theta.size()),
        distribution2(theta.size()),
        dimension    (dimension) {
        assert(dimension > 0);
        // initialize parameters
        dtheta = distribution2.statistics(theta);
        // initialize links
        _links[1] = boost::bind(&categorical_fnode_t::dtheta, this);
}

categorical_fnode_t::categorical_fnode_t(const categorical_fnode_t& categorical_fnode) :
        base_t       (categorical_fnode),
        dtheta       (categorical_fnode.dtheta),
        distribution1(categorical_fnode.distribution1),
        distribution2(categorical_fnode.distribution2),
        dimension    (categorical_fnode.dimension) {
        assert(dimension > 0);
        _links[1] = boost::bind(&categorical_fnode_t::dtheta, this);
}

categorical_fnode_t*
categorical_fnode_t::clone() const {
        return new categorical_fnode_t(*this);
}

bool
categorical_fnode_t::link(const std::string& id, variable_node_i& variable_node) {
        if      (id == "output") return base_t::link(0, variable_node);
        if      (id == "theta")  return base_t::link(1, variable_node);
        else return false;
}

double
categorical_fnode_t::free_energy() const
{
        const vector_t counts    = _links[0]();
        const vector_t log_theta = _links[1]();
        assert(counts.size() == log_theta.size());
        double result = 0.0;

        // statistics & parameters
        for (size_t i = 0; i < counts.size(); i++) {
                result += counts[i]*log_theta[i];
        }

        debug(boost::format("factor node %s:%x computed free energy: %d\n")
              % base_t::name() % this % result);

        return result;
}

bool
categorical_fnode_t::is_conjugate(size_t i, variable_node_i& variable_node) const {
        switch (i) {
        case 0: return variable_node.type() == typeid(categorical_distribution_t);
        case 1: return variable_node.type() == typeid(  dirichlet_distribution_t);
        default: assert(false);
        }
}

const p_message_t&
categorical_fnode_t::operator()(size_t i) {
        switch (i) {
        case 0: return message1();
        case 1: return message2();
        default: assert(false);
        }
}

const p_message_t&
categorical_fnode_t::message1() {
        const vector_t theta = _links[1]();

        debug("categorical message 1 (categorical): " << this->name() << endl);
        distribution1 = categorical_distribution_t(theta);
        debug(endl);

        return distribution1;
}

const p_message_t&
categorical_fnode_t::message2() {
        const vector_t alpha = _links[0]();

        debug("categorical message 2 (dirichlet): " << this->name() << endl);
        distribution2 = dirichlet_distribution_t(alpha);
        debug(endl);

        return distribution1;
}
