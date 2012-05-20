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

// void test_tree1() {
//         cout << "Test 1:" << endl;

//         code_t observations[][2] = {
//                 {0, 0}, {0, 0}, {1, 0}
//         };
//         boost::array<double, alphabet_size> p;
//         p[0] = 0.25;
//         p[1] = 0.25;
//         p[2] = 0.25;
//         p[3] = 0.25;

//         double gradient = 0;
//         for (size_t i = 0; i < 3; i++) {
//                 pt_leaf_t n2( observations[i][0], 0.30, "n2");
//                 pt_leaf_t n3( observations[i][1], 0.30, "n3");
//                 pt_root_t n1(-1, &n2, &n3);

//                 pt_gradient_t<code_t, alphabet_size> result(&n1);
//                 gradient += result.eval(&n2, p);
//         }
//         cout << "gradient for n2: "
//              << gradient
//              << endl;
// }

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

void test_tree3() {

        pt_leaf_t n2( 2, 0.30, "n2");
        pt_leaf_t n3( 2, 0.30, "n3");
        pt_root_t n1(-1, &n2, &n3);

        pt_gradient_t<code_t, alphabet_size> result(&n1);

        cout << result.normalization() << endl;
}

int main(void) {
//        test_tree1();
//        test_tree2();
        test_tree3();

        alignment_t alignment("test.fa");

        return 0.0;
}
