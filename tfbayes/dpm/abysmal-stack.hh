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

#ifndef ABYSMAL_STACK_HH
#define ABYSMAL_STACK_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <iostream>

template <typename T>
class AbysmalStack {
public:
        AbysmalStack(size_t depth)
                : _depth(depth), _bottom(0), _content(depth, 0)
                {}
        
        void push(T letter) {
                _content[_bottom] = letter;
                _bottom = (_bottom + 1) % _depth;
        }

        const T& operator[](size_t pos) const {
                return _content[(_bottom+pos)%_depth];
        }

        T& operator[](size_t pos) {
                return _content[(_bottom+pos)%_depth];
        }

//        friend std::ostream& operator<< <> (std::ostream& o, const AbysmalStack<T>& stack);

protected:
        size_t _depth;
        size_t _bottom;
        std::vector<T> _content;
};

#endif /* ABYSMAL_STACK_HH */
