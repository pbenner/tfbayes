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

#ifndef DPM_GAUSSIAN_INTERFACE_HH
#define DPM_GAUSSIAN_INTERFACE_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <interface.hh>

extern "C" void _dpm_gaussian_init(
        int samples,
        double alpha,
        Bayes::Matrix* _Sigma,
        Bayes::Matrix* _Sigma_0,
        Bayes::Vector* _mu_0,
        Bayes::Vector* _pi);
extern "C" void _dpm_gaussian_sample(unsigned int n, unsigned int burnin);

#endif /* DPM_GAUSSIAN_INTERFACE_HH */
