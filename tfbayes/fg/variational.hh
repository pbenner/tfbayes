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

#ifndef __TFBAYES_FG_VARIATIONAL_HH__
#define __TFBAYES_FG_VARIATIONAL_HH__

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/fg/node-types.hh>
#include <tfbayes/utility/default-operator.hh>

// variable node specializations
////////////////////////////////////////////////////////////////////////////////
class normal_vnode_t : public exponential_vnode_t<normal_distribution_t> {
public:
        typedef exponential_vnode_t<normal_distribution_t> base_t;

        normal_vnode_t(std::string name)
                : base_t(normal_distribution_t(0, 0.1), name)
                { }
};
class gamma_vnode_t : public exponential_vnode_t<gamma_distribution_t> {
public:
        typedef exponential_vnode_t<gamma_distribution_t> base_t;

        gamma_vnode_t(std::string name)
                : base_t(gamma_distribution_t(1, 0.1), name)
                { }
};

class normal_data_t : public data_vnode_t<normal_distribution_t> {
public:
        typedef data_vnode_t<normal_distribution_t> base_t;

        normal_data_t(std::string name)
                : base_t(2, name)
                { }
};

class gamma_data_t : public data_vnode_t<gamma_distribution_t> {
public:
        typedef data_vnode_t<gamma_distribution_t> base_t;

        gamma_data_t(std::string name)
                : base_t(2, name)
                { }
};

// factor node specializations
////////////////////////////////////////////////////////////////////////////////
class normal_fnode_t : public factor_node_t<3> {
public:
        typedef factor_node_t<3> base_t;

        normal_fnode_t(const std::string& name,
                       double mean, double precision, size_t dimension = 1);
        normal_fnode_t(const normal_fnode_t& normal_fnode);

        virtual normal_fnode_t* clone() const;

        friend void swap(normal_fnode_t& left,
                         normal_fnode_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
                swap(left.dmean,         right.dmean);
                swap(left.dprecision,    right.dprecision);
                swap(left.distribution1, right.distribution1);
                swap(left.distribution2, right.distribution2);
                swap(left.distribution3, right.distribution3);
                swap(left.dimension,     right.dimension);
        }
        virtual_assignment_operator(normal_fnode_t);
        derived_assignment_operator(normal_fnode_t, factor_node_i);

        using base_t::link;
        virtual bool link(const std::string& id, variable_node_i& variable_node);
        virtual double free_energy() const;

protected:
        virtual bool is_conjugate(size_t i, variable_node_i& variable_node) const;
        virtual const p_message_t& operator()(size_t i);

        // message preparation
        const p_message_t& message1();
        const p_message_t& message2();
        const p_message_t& message3();

        // parameters
        q_message_t dmean;
        q_message_t dprecision;

        // messages
        normal_distribution_t distribution1;
        normal_distribution_t distribution2;
         gamma_distribution_t distribution3;

        // dimension of the space this node lives on
        size_t dimension;
};

class gamma_fnode_t : public factor_node_t<3> {
public:
        typedef factor_node_t<3> base_t;

        gamma_fnode_t(const std::string& name,
                      double shape, double rate, size_t dimension = 1);
        gamma_fnode_t(const gamma_fnode_t& gamma_fnode);

        virtual gamma_fnode_t* clone() const;

        friend void swap(gamma_fnode_t& left,
                         gamma_fnode_t& right) {
                using std::swap;
                swap(static_cast<base_t&>(left),
                     static_cast<base_t&>(right));
                swap(left.dshape,        right.dshape);
                swap(left.drate,         right.drate);
                swap(left.distribution1, right.distribution1);
                swap(left.distribution3, right.distribution3);
                swap(left.dimension,     right.dimension);
        }
        virtual_assignment_operator(gamma_fnode_t);
        derived_assignment_operator(gamma_fnode_t, factor_node_i);

        using base_t::link;
        virtual bool link(const std::string& id, variable_node_i& variable_node);
        virtual double free_energy() const;

protected:
        virtual bool is_conjugate(size_t i, variable_node_i& variable_node) const;
        virtual const p_message_t& operator()(size_t i);

        // message preparation
        const p_message_t& message1();
        const p_message_t& message3();

        // parameters
        q_message_t dshape;
        q_message_t drate;

        // messages
        gamma_distribution_t distribution1;
        gamma_distribution_t distribution3;

        // dimension of the space this node lives on
        size_t dimension;
};

#endif /* __TFBAYES_FG_VARIATIONAL_HH__ */
