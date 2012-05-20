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

class pmut_t {
public:
        pmut_t(pt_node_t* node, bool mutation = true)
                : node(node), mutation(mutation) { }

        operator bool() const {
                return mutation;
        }
        bool operator==(const pmut_t& m) const {
                return node == m.node && mutation == m.mutation;
        }
        bool operator<(const pmut_t& m) const {
                if (node < m.node) {
                        return true;
                }
                if (node == m.node) {
                        return mutation < m.mutation;
                }
                return false;
        }

        pt_node_t* node;

private:
        bool mutation;
};

size_t
hash_value(const pmut_t& pmut)
{
        return (size_t)pmut.node;
}

#include <vector>
#include <set>

class mutation_product_t : public std::multiset<pmut_t> {
public:
        mutation_product_t()
                : std::multiset<pmut_t>() {}
        mutation_product_t(const pmut_t& mutation)
                : std::multiset<pmut_t>() {
                insert(mutation);
        }

        mutation_product_t operator*=(const pmut_t& mutation) {
                insert(mutation);
                return *this;
        }
        mutation_product_t operator*=(const mutation_product_t& product) {

                for (mutation_product_t::const_iterator it = product.begin(); it != product.end(); it++) {
                        operator*=(*it);
                }
                return *this;
        }
};

class mutation_coefficient_t : public boost::unordered_map<mutation_product_t, double> {
public:
        mutation_coefficient_t()
                : boost::unordered_map<mutation_product_t, double>() {}
        mutation_coefficient_t(const pmut_t& mutation)
                : boost::unordered_map<mutation_product_t, double>() {
                operator[](mutation) += 1.0;
        }
        mutation_coefficient_t(const mutation_product_t& product)
                : boost::unordered_map<mutation_product_t, double>() {
                operator[](product) += 1.0;
        }

        operator bool() const {
                return size();
        }
        mutation_coefficient_t operator+=(
                const mutation_coefficient_t& coefficient) {

                for (mutation_coefficient_t::const_iterator it = coefficient.begin(); it != coefficient.end(); it++) {
                        operator[](it->first) += it->second;
                        if (operator[](it->first) == 0.0) {
                                erase(it->first);
                        }
                }
                return *this;
        }
        mutation_coefficient_t operator*=(
                const mutation_coefficient_t& coefficient) {

                mutation_coefficient_t tmp;
                for (mutation_coefficient_t::const_iterator it = begin(); it != end(); it++) {
                        for (mutation_coefficient_t::const_iterator is = coefficient.begin(); is != coefficient.end(); is++) {
                                mutation_product_t product(it->first);
                                product *= (is->first);
                                tmp[product] += (it->second)*(is->second);
                        }
                }
                operator=(tmp);

                return *this;
        }
};

ostream& operator<< (ostream& o, const mutation_product_t product) {
        for (mutation_product_t::const_iterator it = product.begin(); it != product.end(); it++) {
                const pmut_t& mutation = *it;
                if (mutation) {
                        o << "(1-M(" << mutation.node->name << ")) ";
                }
                else {
                        o << "M(" << mutation.node->name << ") ";
                }
        }
        return o;
}

ostream& operator<< (ostream& o, const mutation_coefficient_t& coefficient) {

        o << "[ ";
        for (mutation_coefficient_t::const_iterator it = coefficient.begin(); it != coefficient.end(); it++) {
                if (it != coefficient.begin()) {
                        o << "+ ";
                }
                if (it->second != 1.0) {
                        o << it->second << " ";
                }
                o << it->first;
        }
        o << "]";
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


        mutation_product1 *= pmut_t(&n2, true);
        mutation_product1 *= pmut_t(&n3, false);

        mutation_product2 *= pmut_t(&n2, false);
        mutation_product2 *= pmut_t(&n3, false);

        mutation_coefficient += mutation_product1;
        mutation_coefficient += mutation_product2;

        polynomial_term_t<code_t, alphabet_size, mutation_coefficient_t> term1(mutation_coefficient);
        term1.exponent()[2] += 2;
        polynomial_t     <code_t, alphabet_size, mutation_coefficient_t> poly1(term1);

        mutation_coefficient *= pmut_t(&n1, true);

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
