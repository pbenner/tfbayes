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

#include <node-types.hh>

// variable node specializations
////////////////////////////////////////////////////////////////////////////////
typedef variable_node_t<normal_distribution_t> normal_vnode_t;
typedef variable_node_t< gamma_distribution_t>  gamma_vnode_t;

// factor node specializations
////////////////////////////////////////////////////////////////////////////////
class normal_fnode_t : public factor_node_t<3> {
public:
        typedef factor_node_t<3> base_t;

        normal_fnode_t(double mean, double precision);

        virtual normal_fnode_t* clone() const;

        virtual bool link(const std::string& id, variable_node_i& variable_node);
protected:
        virtual bool is_conjugate(size_t i, variable_node_i& variable_node) const;
        virtual const p_message_t& initial_message(size_t i) const;
        virtual const p_message_t& message(size_t i);

        // message preparation
        const p_message_t& message1();
        const p_message_t& message2();
        const p_message_t& message3();

        // parameters
        dirac_distribution_t dmean;
        dirac_distribution_t dprecision;

        // messages
        normal_distribution_t distribution1;
        normal_distribution_t distribution2;
         gamma_distribution_t distribution3;
};

class gamma_fnode_t : public factor_node_t<3> {
public:
        typedef factor_node_t<3> base_t;

        gamma_fnode_t(double shape, double rate);

        virtual gamma_fnode_t* clone() const;

        virtual bool link(const std::string& id, variable_node_i& variable_node);
protected:
        virtual bool is_conjugate(size_t i, variable_node_i& variable_node) const;
        virtual const p_message_t& initial_message(size_t i) const;
        virtual const p_message_t& message(size_t i);

        // message preparation
        const p_message_t& message1();
        const p_message_t& message3();

        // parameters
        dirac_distribution_t dshape;
        dirac_distribution_t drate;

        // messages
        gamma_distribution_t distribution1;
        gamma_distribution_t distribution3;
};

#endif /* __TFBAYES_FG_VARIATIONAL_HH__ */
