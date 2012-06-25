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
#include <cassert>

#include <sys/time.h>

#include <phylotree-polynomial.hh>
#include <utility.hh>

using namespace std;

#define alphabet_size 4
typedef short code_t;

class pt_test_t : public pt_root_t {
public:
        pt_test_t(size_t n)
                : pt_root_t(-1) {

                assert(n >= 1);

                pt_node_t* pt_last = this;

                for (size_t i = 1; i < n; i++) {
                        pt_last->left  = new pt_node_t(-1, 1.0);
                        pt_last->right = new pt_node_t(-1, 1.0);
                        pt_last = pt_last->left;
                }
                pt_node_t::id_t i = set_id(this, 0)+1;
                node_map = node_map_t(i, (pt_node_t*)NULL);
                create_map(this);
        }
};

void random_init(pt_node_t* pt_node)
{
        if (pt_node->leaf()) {
                pt_node->x = rand()%alphabet_size;
        }
        else {
                random_init(pt_node->left);
                random_init(pt_node->right);
        }
}

double expected_terms(size_t n)
{
        pt_root_t* pt_root = new pt_test_t(n);
        double result = 0;

        for (size_t i = 0; i < 1000; i++)
        {
                random_init(pt_root);
                pt_simple_polynomial_t<code_t, alphabet_size> poly(pt_root);
//                pt_polynomial_t<code_t, alphabet_size> poly(pt_root);
                result = (i*(double)result + (double)poly.size())/((double)i+1.0);

        }
        pt_root->destroy();

        return result;
}

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

int main(void)
{
        init();

//        size_t N = 20;
        size_t N = 6;

        cout << "c(";
        for (size_t n = 1; n <= N; n++) {
                cout << 2*n-1;
                if (n < N) {
                     cout << ", ";
                }
        }
        cout << ")" << endl;
        cout << "c(";
        for (size_t n = 1; n <= N; n++) {
                double result = expected_terms(n);
                cout << result;
                if (n < N) {
                     cout << ", ";
                }
        }
        cout << ")" << endl;

        return 0.0;
}
