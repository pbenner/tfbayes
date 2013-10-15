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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <tfbayes/phylotree/utility.hh>

using namespace std;

ostream& operator<< (ostream& o, const exponent_t<4, char>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<4, char>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<4, char>& polynomial) {
        for (polynomial_t<4, char>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}
ostream& operator<< (ostream& o, const exponent_t<5, char>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];
        if(exponent[4]) o << " Pn^" << exponent[4];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<5, char>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<5, char>& polynomial) {
        for (polynomial_t<5, char>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}

////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const exponent_t<4, short>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<4, short>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<4, short>& polynomial) {
        for (polynomial_t<4, short>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}
ostream& operator<< (ostream& o, const exponent_t<5, short>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];
        if(exponent[4]) o << " Pn^" << exponent[4];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<5, short>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<5, short>& polynomial) {
        for (polynomial_t<5, short>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}

////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const exponent_t<4, float>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<4, float>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<4, float>& polynomial) {
        for (polynomial_t<4, float>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}
ostream& operator<< (ostream& o, const exponent_t<5, float>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];
        if(exponent[4]) o << " Pn^" << exponent[4];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<5, float>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<5, float>& polynomial) {
        for (polynomial_t<5, float>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}

////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& o, const exponent_t<4, double>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<4, double>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<4, double>& polynomial) {
        for (polynomial_t<4, double>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}
ostream& operator<< (ostream& o, const exponent_t<5, double>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];
        if(exponent[4]) o << " Pn^" << exponent[4];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<5, double>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<5, double>& polynomial) {
        for (polynomial_t<5, double>::const_iterator it = polynomial.begin();
             it != polynomial.end(); it++) {
                if (it != polynomial.begin()) {
                        o << " + " << *it;
                }
                else {
                        o << *it;
                }
        }

        return o;
}

////////////////////////////////////////////////////////////////////////////////

size_t hash_value(const exponent_t<2, char>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<3, char>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<4, char>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}


size_t hash_value(const exponent_t<5, char>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<10, char>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}

size_t hash_value(const exponent_t<2, short>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<3, short>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<4, short>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}


size_t hash_value(const exponent_t<5, short>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<10, short>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}

size_t hash_value(const exponent_t<2, float>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<3, float>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<4, float>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}

size_t hash_value(const exponent_t<5, float>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<10, float>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}

size_t hash_value(const exponent_t<2, double>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;

        return seed;
}

size_t hash_value(const exponent_t<3, double>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;

        return seed;
}

size_t hash_value(const exponent_t<4, double>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;

        return seed;
}

size_t hash_value(const exponent_t<5, double>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;

        return seed;
}

size_t hash_value(const exponent_t<10, double>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] <<  0;
        seed += (size_t)exponent[1] <<  2;
        seed += (size_t)exponent[2] <<  4;
        seed += (size_t)exponent[3] <<  6;
        seed += (size_t)exponent[4] <<  8;
        seed += (size_t)exponent[5] << 10;
        seed += (size_t)exponent[6] << 12;
        seed += (size_t)exponent[7] << 14;
        seed += (size_t)exponent[8] << 16;
        seed += (size_t)exponent[9] << 18;

        return seed;
}
