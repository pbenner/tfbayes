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

#ifndef __TFBAYES_UTILITY__PROGRESS_HH___
#define __TFBAYES_UTILITY__PROGRESS_HH___

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <sstream>
#include <iomanip>

class progress_t : public std::string {
public:
        progress_t(double p) : std::string() {
                std::stringstream ss;

                // carriage return
                ss << "\33[2K\r" << "|";
                for (size_t i = 1; i < LINE_WIDTH-1; i++) {
                        if (i/(double)LINE_WIDTH < p) {
                                ss << ">";
                        }
                        else {
                                ss << " ";
                        }
                }
                ss << "| " << std::fixed << std::setprecision(2)
                   << p*100 << "%";
                // add newline if finished
                if (p == 1.0) {
                        ss << "\n";
                }
                // copy to string
                std::string::operator=(ss.str());
        }
protected:
        static const ssize_t LINE_WIDTH = 40;
};

#endif /* __TFBAYES_UTILITY__PROGRESS_HH___ */
