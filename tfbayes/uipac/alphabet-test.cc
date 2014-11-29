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


#include <iostream>

#include <tfbayes/uipac/alphabet.hh>

using namespace std;

void
test_nucleotide_code()
{
        alphabet_t b = nucleotide_alphabet_t();
        alphabet_t a = b;

        for (alphabet_code_t i = 0; i < 126; i++) {
                cout << i << " is coded as " << (int)a.code(i) << " and decoded as " << a.decode(a.code(i));
                if (a.element(i)) {
                        cout << " (is element)";
                }
                cout << endl;
        }
}

void
test_protein_code()
{
        alphabet_t b = protein_alphabet_t();
        alphabet_t a = b;

        for (alphabet_code_t i = 0; i < 126; i++) {
                cout << i << " is coded as " << (int)a.code(i) << " and decoded as " << a.decode(a.code(i));
                if (a.element(i)) {
                        cout << " (is element)";
                }
                cout << endl;
        }
}

int
main(void)
{
        test_nucleotide_code();
        cout << endl;
        test_protein_code();

        return 0;
}
