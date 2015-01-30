/* Copyright (C) 2015 Philipp Benner
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
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/utility/probability.hh>
#include <tfbayes/entropy/entropy-distribution.hh>

std::ostream& operator<<(std::ostream& o, const std::vector<probability_t>& v) {
        for (size_t i = 0; i < v.size(); i++) {
                if (i != 0) o << " ";
                o << std::fixed << std::setprecision(8) << v[i];
        }
        return o;
}

template <typename T>
void seed_rng(T& rng)
{
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

using namespace std;

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s SAMPLES A1 A2\n\n", pname);
}

int
main(int argc, char *argv[])
{
        if (argc != 4) {
                print_usage(argv[0], stderr);
                exit(EXIT_FAILURE);
        }
        size_t n  = atoi(argv[1]);
        double a1 = atof(argv[2]);
        double a2 = atof(argv[3]);

        boost::random::mt19937 gen; seed_rng(gen);
        boost::random::entropy_distribution<double, probability_t> dist(3, a1, a2);

        for (size_t i = 0; i < n; i++) {
                cout << dist(gen, 0.1)
                     << endl;
        }
        cerr << "acceptance ratio: "
             << dist.acceptance_ratio()
             << endl;
}
