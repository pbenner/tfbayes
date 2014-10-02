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

#ifndef __TFBAYES_UIPAC_ALPHABET_HH__
#define __TFBAYES_UIPAC_ALPHABET_HH__

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
        alphabet_code_t complement(alphabet_code_t code) const {
                return _complement[static_cast<size_t>(code)];
        }
        bool element(alphabet_code_t letter) const {
                return _code[static_cast<size_t>(letter)] >= 0;
        }
        size_t size() const {
                return _size;
        }
        bool operator==(const alphabet_t& alphabet) const {
                return _code == alphabet._code;
        }
        bool operator!=(const alphabet_t& alphabet) const {
                return _code != alphabet._code;
        }
protected:
        // a derived class may use this constructor to
        // initialize the codebooks
        alphabet_t(const alphabet_code_t* code,
                   const alphabet_code_t* decode,
                   const alphabet_code_t* complement,
                   size_t size)
                : _code(code), _decode(decode), _complement(complement), _size(size)
                { }
        // codebook pointer
        const alphabet_code_t*       _code;
        const alphabet_code_t*     _decode;
        const alphabet_code_t* _complement;
        // the size of the alphabet
        size_t _size;
};

class nucleotide_alphabet_t : public alphabet_t {
public:
        nucleotide_alphabet_t() 
                : alphabet_t(_code, _decode, _complement, 5)
                { }
protected:
        // codebooks specific to this class
        static const alphabet_code_t       _code[];
        static const alphabet_code_t     _decode[];
        static const alphabet_code_t _complement[];
};

#endif /* __TFBAYES_UIPAC_ALPHABET_HH__ */
