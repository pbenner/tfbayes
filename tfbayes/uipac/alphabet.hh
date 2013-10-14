/* Copyright (C) 2011-2013 Philipp Benner
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

#ifndef ALPHABET_HH
#define ALPHABET_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cstddef>

typedef char alphabet_code_t;

class alphabet_t {
public:
        alphabet_code_t code(alphabet_code_t letter) const {
                return _code[static_cast<size_t>(letter)];
        }
        alphabet_code_t decode(alphabet_code_t code) const {
                return _decode[static_cast<size_t>(code)];
        }
        bool element(alphabet_code_t letter) const {
                return _code[static_cast<size_t>(letter)] >= 0;
        }
        size_t size() const {
                return _size;
        }
protected:
        // a derived class may use this constructor to
        // initialize the codebooks
        alphabet_t(const alphabet_code_t* code,
                   const alphabet_code_t* decode,
                   size_t size)
                : _code(code), _decode(decode), _size(size)
                { }
        // codebook pointer
        const alphabet_code_t*   _code;
        const alphabet_code_t* _decode;
        // the size of the alphabet
        size_t _size;
};

class nucleotide_alphabet_t : public alphabet_t {
public:
        nucleotide_alphabet_t() 
                : alphabet_t(_code, _decode, 5)
                { }
protected:
        // codebooks specific to this class
        static const alphabet_code_t   _code[];
        static const alphabet_code_t _decode[];
};

#endif /* ALPHABET_HH */
