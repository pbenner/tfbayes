/* Copyright (C) 2012 Philipp Benner
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

#ifndef STRTOOLS_HH
#define STRTOOLS_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

static __inline__
const std::string strip(const std::string& str) {
        const std::string& whitespace = " \t\n";
        const size_t begin = str.find_first_not_of(whitespace);
        if (begin == std::string::npos) {
                        return "";
        }
        const size_t end   = str.find_last_not_of(whitespace);
        const size_t range = end - begin + 1;
        
        return str.substr(begin, range);
}

static __inline__
std::vector<std::string> token(const std::string& str, char t) {
        std::string token;
        std::vector<std::string> tokens;
        std::istringstream iss(str);
        while (getline(iss, token, t)) {
                tokens.push_back(strip(token));
        }
        return tokens;
}

/* introduce a new line every n characters */
static __inline__
std::string split_string(const std::string& str, size_t n)
{
        std::string result(str);
        for (std::string::iterator it = result.begin()+n;
             it < result.end(); it+=n) {
                it = result.insert(it, '\n');
                it++;
        }
        return result;
}

#endif /* STRTOOLS_HH */
