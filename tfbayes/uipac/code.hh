/* Copyright (C) 2011, 2012 Philipp Benner
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

#ifndef CODE_HH
#define CODE_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

bool is_nucleotide(const char S);

template <typename CODE_TYPE>
CODE_TYPE code_nucleotide(char a)
{
        switch (a) {
        case 'A':
        case 'a':
                return 1;
        case 'C':
        case 'c':
                return 3;
        case 'G':
        case 'g':
                return 0;
        case 'T':
        case 't':
                return 2;
        case '-':
                return 4;
        case 'N':
        case 'n':
                return -1;
        default:
                break;
        }
        return -1;
}

template <typename CODE_TYPE>
char decode_nucleotide(CODE_TYPE a)
{
        switch (a) {
        case 1:
                return 'A';
        case 3:
                return 'C';
        case 0:
                return 'G';
        case 2:
                return 'T';
        case 4:
                return '-';
        case -1:
                return 'N';
        default:
                break;
        }
        return 'X';
}

#endif /* CODE_HH */
