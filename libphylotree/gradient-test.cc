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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>

#include <phylotree-gradient.hh>
#include <utility.hh>

using namespace std;

void test_tree1() {
        cout << "Test 1:" << endl;

        code_t observations[][2] = {
                {0, 0}, {0, 0}, {1, 0}
        };
        boost::array<double, alphabet_size> p;
        p[0] = 0.25;
        p[1] = 0.25;
        p[2] = 0.25;
        p[3] = 0.25;

        double gradient = 0;
        for (size_t i = 0; i < 3; i++) {
                pt_leaf_t n2( observations[i][0], 0.30, "n2");
                pt_leaf_t n3( observations[i][1], 0.30, "n3");
                pt_root_t n1(-1, &n2, &n3);

                pt_gradient_t<code_t, alphabet_size> result(&n1);
                gradient += result.eval(&n2, p);
        }
        cout << "gradient for n2: "
             << gradient
             << endl;
}

#include <tfbayes/fasta.hh>

static
const string strip(const std::string& str)
{
        const std::string& whitespace = " \t\n";
        const size_t begin = str.find_first_not_of(whitespace);
        if (begin == string::npos) {
                return "";
        }
        const size_t end   = str.find_last_not_of(whitespace);
        const size_t range = end - begin + 1;

        return str.substr(begin, range);
}

static
vector<string> token(const string& str, char t) {
        string token;
        vector<string> tokens;
        istringstream iss(str);
        while (getline(iss, token, t)) {
                tokens.push_back(strip(token));
        }
        return tokens;
}

#include <set>
#include <boost/unordered_map.hpp>

class alignment_t : public boost::unordered_map<string, string> {
public:
        alignment_t(const char* filename)
                : boost::unordered_map<string, string>() {

                FastaParser parser(filename);
                string sequence;

                while ((sequence = parser.read_sequence()) != "") {
                        string taxon = token(parser.description()[0], '.')[0];
                        taxa.insert(taxon);
                        operator[](taxon) = sequence;
                }
        }

        set<string> taxa;
};

class mutation_t {
public:
        mutation_t(pt_node_t* node, bool mutation = true)
                : node(node), mutation(mutation) { }

        operator bool() const {
                return mutation;
        }

        pt_node_t* node;

private:
        bool mutation;
};

#include <vector>

class mutation_product_t : public std::vector<mutation_t> {
public:
        mutation_product_t()
                : std::vector<mutation_t>() {}
        mutation_product_t(const mutation_t& mutation)
                : std::vector<mutation_t>() {
                push_back(mutation);
        }

        mutation_product_t operator*=(const mutation_t& mutation) {
                push_back(mutation);
                return *this;
        }
        mutation_product_t operator*=(const mutation_product_t& product) {

                for (mutation_product_t::const_iterator it = product.begin(); it != product.end(); it++) {
                        operator*=(*it);
                }
                return *this;
        }
};

class mutation_coefficient_t : public std::vector<mutation_product_t> {
public:
        mutation_coefficient_t()
                : std::vector<mutation_product_t>() {}
        mutation_coefficient_t(const mutation_t& mutation)
                : std::vector<mutation_product_t>() {
                push_back(mutation);
        }
        mutation_coefficient_t(const mutation_product_t& product)
                : std::vector<mutation_product_t>() {
                push_back(product);
        }

        operator bool() const {
                return size();
        }
        mutation_coefficient_t operator+=(
                const mutation_product_t& product) {
                push_back(product);
                return *this;
        }
        mutation_coefficient_t operator+=(
                const mutation_coefficient_t& coefficient) {

                for (mutation_coefficient_t::const_iterator it = coefficient.begin(); it != coefficient.end(); it++) {
                        operator+=(*it);
                }
                return *this;
        }
        mutation_coefficient_t operator*=(
                const mutation_coefficient_t& coefficient) {

                for (mutation_coefficient_t::iterator it = begin(); it != end(); it++) {
                        for (mutation_coefficient_t::const_iterator is = coefficient.begin(); is != coefficient.end(); is++) {
                                (*it) *= (*is);
                        }
                }
                return *this;
        }
};

ostream& operator<< (ostream& o, const mutation_product_t product) {
        for (mutation_product_t::const_iterator it = product.begin(); it != product.end(); it++) {
                const mutation_t& mutation = *it;
                if (mutation) {
                        o << "(1-e^" << mutation.node->name << ") ";
                }
                else {
                        o << "e^" << mutation.node->name << " ";
                }
        }
        return o;
}

ostream& operator<< (ostream& o, const mutation_coefficient_t coefficient) {

        o << "[ ";
        for (mutation_coefficient_t::const_iterator it = coefficient.begin(); it != coefficient.end(); it++) {
                if (it != coefficient.begin()) {
                        o << "+ ";
                }
                o << *it;
        }
        o << "]";
        return o;
}

ostream& operator<< (ostream& o, const exponent_t<code_t, alphabet_size>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<code_t, alphabet_size, mutation_coefficient_t>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<code_t, alphabet_size, mutation_coefficient_t>& polynomial) {
        for (polynomial_t<code_t, alphabet_size, mutation_coefficient_t>::const_iterator it = polynomial.begin();
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


void test_tree2() {
        cout << "Test 2:" << endl;

        pt_leaf_t n2( 1, 0.20, "n2");
        pt_leaf_t n3( 2, 0.30, "n3");
        pt_root_t n1(-1, &n2, &n3, "n1");

        mutation_product_t mutation_product1;
        mutation_product_t mutation_product2;
        mutation_coefficient_t mutation_coefficient;

        mutation_product1.push_back(mutation_t(&n2, true));
        mutation_product1.push_back(mutation_t(&n3, false));

        mutation_product2.push_back(mutation_t(&n2, false));
        mutation_product2.push_back(mutation_t(&n3, false));

        mutation_coefficient += mutation_product1;
        mutation_coefficient += mutation_product2;

        polynomial_term_t<code_t, alphabet_size, mutation_coefficient_t> term1(mutation_coefficient);
        term1.exponent()[2] += 2;
        polynomial_t     <code_t, alphabet_size, mutation_coefficient_t> poly1(term1);

        mutation_coefficient *= mutation_t(&n1, true);

        polynomial_term_t<code_t, alphabet_size, mutation_coefficient_t> term2(mutation_coefficient);
        term2.exponent()[2] += 2;
        polynomial_t     <code_t, alphabet_size, mutation_coefficient_t> poly2(term2);

        cout << poly1 << endl;
        cout << poly2 << endl;
        cout << poly1*poly2 << endl;
}

int main(void) {
        test_tree1();
        test_tree2();

        alignment_t alignment("test.fa");


        return 0.0;
}
