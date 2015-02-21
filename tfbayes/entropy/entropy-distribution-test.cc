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

#include <getopt.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <tfbayes/utility/probability.hh>
#include <tfbayes/entropy/entropy-distribution.hh>
#include <tfbayes/utility/strtools.hh>

typedef probability_t<> p_t;
typedef std::vector<p_t> p_vector_t;

std::ostream& operator<<(std::ostream& o, const p_vector_t& v) {
        for (size_t i = 0; i < v.size(); i++) {
                if (i != 0) o << " ";
                o << std::fixed << std::setprecision(8) << v[i];
        }
        return o;
}

using namespace std;

// options
////////////////////////////////////////////////////////////////////////////////

struct options_t {
        double a1;
        double a2;
        bool print_entropy;
        bool verbose;
        // default constructor
        options_t()
                : a1            (1.0)
                , a2            (1.0)
                , print_entropy (false)
                , verbose       (false)
                { }
};

static options_t options;

// usage and version information
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] SAMPLES CARDINALITY\n\n", pname);
        (void)fprintf(fp,
                      "Options:\n"
                      "      --alpha=float             - parameter of the dirichlet proposal distribution\n"
                      "      --counts=A1:A2            - pseudocounts for the beta prior\n"
                      "      --print-entropy           - print entropy instead of the distributions\n"
                      "      --variance=float          - variance of the proposal distribution\n"
                      "   -v                           - be verbose\n"
                      "\n"
                      "      --help                    - print help and exit\n"
                      "      --version                 - print version information and exit\n\n");
}

static
void wrong_usage(char *pname, const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `%s --help' for more information.\n", pname);

        exit(EXIT_FAILURE);

}

static
void print_version(FILE *fp)
{
        (void)fprintf(fp,
                      "This is free software, and you are welcome to redistribute it\n"
                      "under certain conditions; see the source for copying conditions.\n"
                      "There is NO warranty; not even for MERCHANTABILITY or FITNESS\n"
                      "FOR A PARTICULAR PURPOSE.\n\n");
}

// main program
////////////////////////////////////////////////////////////////////////////////

template <typename T>
void seed_rng(T& rng)
{
        struct timeval tv;
        gettimeofday(&tv, NULL);
        rng.seed(tv.tv_sec*tv.tv_usec);
}

void
sample(size_t samples, size_t k)
{
        boost::random::mt19937 gen; seed_rng(gen);
        entropy_distribution_t<double, p_t> d(k, options.a1, options.a2);

        if (options.print_entropy) {
                for (size_t i = 0; i < samples; i++) {
                        cout << entropy(d(gen))
                             << endl;
                }
        }
        else {
                for (size_t i = 0; i < samples; i++) {
                        cout << d(gen)
                             << endl;
                }
        }
        if (options.verbose) {
                cerr << "acceptance ratio: "
                     << d.acceptance_ratio()
                     << endl;
        }
}

// main function
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
        extern int optind;
        vector<string> tokens;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "counts",               1, 0, 'b' },
                        { "print-entropy",        0, 0, 'c' },
                        { "help",                 0, 0, 'h' },
                        { "version",              0, 0, 'x' },
                        { 0,                      0, 0,  0  }
                };

                c = getopt_long(argc, argv, "v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'b':
                        tokens = token(string(optarg), ':');
                        if (tokens.size() != 2) {
                                wrong_usage(argv[0], "Invalid number of pseudocounts.");
                        }
                        options.a1 = atof(tokens[0].c_str());
                        options.a2 = atof(tokens[1].c_str());
                        break;
                case 'c':
                        options.print_entropy = true;
                        break;
                case 'v':
                        options.verbose = true;
                        break;
                case 'h':
                        print_usage(argv[0], stdout);
                        exit(EXIT_SUCCESS);
                case 'x':
                        print_version(stdout);
                        exit(EXIT_SUCCESS);
                default:
                        wrong_usage(argv[0], NULL);
                        exit(EXIT_FAILURE);
                }
        }
        if(optind+2 != argc) {
                wrong_usage(argv[0], "Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }
        sample(atoi(argv[optind]), atoi(argv[optind+1]));

        return 0;
}
