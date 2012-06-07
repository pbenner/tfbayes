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

#include <alignment.hh>
#include <phylotree.hh>
#include <phylotree-parser.hh>
#include <phylotree-polynomial.hh>
#include <phylotree-sampler.hh>
#include <phylotree-gradient.hh>
#include <marginal-likelihood.hh>

#include <tfbayes/exception.h>

#include <getopt.h>

#define alphabet_size 5
typedef float code_t;

using namespace std;

ostream& operator<< (ostream& o, const exponent_t<code_t, alphabet_size>& exponent) {
        if(exponent[0]) o << " Pa^" << exponent[0];
        if(exponent[1]) o << " Pc^" << exponent[1];
        if(exponent[2]) o << " Pg^" << exponent[2];
        if(exponent[3]) o << " Pt^" << exponent[3];

        return o;
}
ostream& operator<< (ostream& o, const polynomial_term_t<code_t, alphabet_size>& term) {
        o << term.coefficient()
          << term.exponent();

        return o;
}
ostream& operator<< (ostream& o, const polynomial_t<code_t, alphabet_size>& polynomial) {
        for (polynomial_t<code_t, alphabet_size>::const_iterator it = polynomial.begin();
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

size_t hash_value(const exponent_t<code_t, alphabet_size>& exponent) {
        size_t seed = 0;
        seed += (size_t)exponent[0] << 0;
        seed += (size_t)exponent[1] << 2;
        seed += (size_t)exponent[2] << 4;
        seed += (size_t)exponent[3] << 6;
        seed += (size_t)exponent[4] << 8;

        return seed;
}

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        double alpha;
        size_t burnin;
        double lambda;
        double r;
        size_t max_steps;
        double min_change;
        double epsilon;
        double sigma;
        _options_t()
                : alpha(0.2),
                  burnin(1000),
                  lambda(0.1),
                  r(2.0),
                  max_steps(1000),
                  min_change(0.0005),
                  epsilon(0.001),
                  sigma(0.01)
                { }
} options_t;

static options_t options;

// Main
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] METHOD TREE ALIGNMENT\n\n", pname);
        (void)fprintf(fp,
                      "Methods: gradient-ascent, metropolis-hastings\n"
                      "\n"
                      "Options:\n"
                      "             -a DOUBLE       - pseudo count\n"
                      "             -b INTEGER      - number of burn-in samples for metropolis-hastings\n"
                      "             -m INTEGER      - maximum number of steps\n"
                      "             -n DOUBLE       - stop gradient ascent if change is smaller than this value\n"
                      "             -e DOUBLE       - gradient ascent step size\n"
                      "             -r DOUBLE       - r parameter for the gamma distribution\n"
                      "             -l DOUBLE       - lambda parameter for the gamma distribution\n"
                      "             -s DOUBLE       - sigma^2 parameter for proposal distribution\n"
                      "\n"
                      "   --help                    - print help and exit\n"
                      "   --version                 - print version information and exit\n\n");
}

void wrong_usage(const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `optimization-test --help' for more information.\n");

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

extern FILE *yyin;

pt_root_t* parse_tree_file(const char* file_tree)
{
        yyin = fopen(file_tree, "r");
        if (yyin == NULL) {
                std_err(PERR, "Could not open phylogenetic tree");
        }
        yyparse();
        fclose(yyin);

        return (pt_root_t*)pt_parsetree->convert();
}

void run_optimization(const string& method, const char* file_tree, const char* file_alignment)
{
        /* phylogenetic tree */
        pt_root_t* pt_root = parse_tree_file(file_tree);
        pt_root->init(alphabet_size);

        /* pseudo counts */
        exponent_t<code_t, alphabet_size> alpha;
        for (size_t i = 0; i < alphabet_size; i++) {
                alpha[i] = options.alpha;
        }

        /* alignment */
        alignment_t<code_t> alignment(file_alignment, pt_root);

        if (method == "gradient-ascent") {
                pt_gradient_ascent_t<code_t, alphabet_size> pt_gradient_ascent(pt_root, alignment, alpha, options.r, options.lambda, options.epsilon);
                pt_gradient_ascent.run(options.max_steps, options.min_change);
        }
        else if (method == "metropolis-hastings") {
                normal_jump_t jump(options.sigma);
//                gamma_jump_t jump(1.6, 0.4);
                pt_metropolis_hastings_t<code_t, alphabet_size> pt_metropolis_hastings(pt_root, alignment, alpha, options.r, options.lambda, jump, 0.5);
                pt_metropolis_hastings.burnin(options.burnin);
                pt_metropolis_hastings.sample(options.max_steps);
                //pt_pmcmc_hastings_t<code_t, alphabet_size> pmcmc(10, pt_metropolis_hastings);
                //pmcmc.sample(options.max_steps, options.burnin);
        }
        else {
                cerr << "Unknown optimization method: " << method
                     << endl;
                exit(EXIT_FAILURE);
        }

        stringstream ss;
        pt_root->print(ss, false);
        cout << ss.str()
             << endl;
}

int main(int argc, char *argv[])
{
        const char* file_tree;
        const char* file_alignment;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'v' }
                };

                c = getopt_long(argc, argv, "a:b:e:m:n:r:l:s:hv",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        options.alpha = atof(optarg);
                        break;
                case 'b':
                        options.burnin = atoi(optarg);
                        break;
                case 'e':
                        options.epsilon = atof(optarg);
                        break;
                case 'm':
                        options.max_steps = atoi(optarg);
                        break;
                case 'n':
                        options.min_change = atof(optarg);
                        break;
                case 'r':
                        options.r = atof(optarg);
                        break;
                case 'l':
                        options.lambda = atof(optarg);
                        break;
                case 's':
                        options.sigma = atof(optarg);
                        break;
                case 'h':
                        print_usage(argv[0], stdout);
                        exit(EXIT_SUCCESS);
                case 'v':
                        print_version(stdout);
                        exit(EXIT_SUCCESS);
                default:
                        wrong_usage(NULL);
                        exit(EXIT_FAILURE);
                }
        }
        if(optind+3 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }

        string method(argv[optind]);
        file_tree      = argv[optind+1];
        file_alignment = argv[optind+2];

        run_optimization(method, file_tree, file_alignment);

        return 0;
}
