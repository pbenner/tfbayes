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
typedef variable_node_t<normal_distribution_t> normal_vnode_t;
typedef variable_node_t< gamma_distribution_t>  gamma_vnode_t;

// a factor node is an exponential family with the ability to send
// and receive messages
class normal_fnode_t : public factor_node_t<3> {
public:
        typedef factor_node_t<3> base_t;

        virtual normal_fnode_t* clone() const {
                return new normal_fnode_t(*this);
        }

protected:
        void init(size_t i) {
                switch (i) {
                case '0': outbox[i]->replace(&distribution1);
                case '1': outbox[i]->replace(&distribution2);
                case '2': outbox[i]->replace(&distribution3);
                default: assert(false);
                }
        }
        const p_message_t& message(size_t i) {
                switch (i) {
                case '0': return message1();
                case '1': return message2();
                case '2': return message3();
                default: assert(false);
                }
        }
        const p_message_t& message1() {
                return distribution1;
        }
        const p_message_t& message2() {
                return distribution2;
        }
        const p_message_t& message3() {
                return distribution3;
        }

        // messages
        normal_distribution_t distribution1;
        normal_distribution_t distribution2;
         gamma_distribution_t distribution3;

        normal_distribution_t msg1;
};

#endif /* __TFBAYES_FG_VARIATIONAL_HH__ */
