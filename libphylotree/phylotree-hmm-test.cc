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
#include <iomanip>
#include <string>

#include <phylotree-hmm.hh>
#include <tfbayes/exception.h>
#include <tfbayes/strtools.hh>

#include <getopt.h>

#define alphabet_size 5
typedef float code_t;

using namespace std;

// Input/Output
////////////////////////////////////////////////////////////////////////////////

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
        size_t dimension;
        phylotree_hmm_t<code_t, alphabet_size>::priors_t priors;
        phylotree_hmm_t<code_t, alphabet_size>::matrix_t transition;
        bool verbose;
        _options_t()
                : dimension(2),
                  priors(),
                  transition(),
                  verbose(false) {

                exponent_t<code_t, alphabet_size> alpha_0;
                exponent_t<code_t, alphabet_size> alpha_1;
                for (size_t i = 0; i < alphabet_size; i++) {
                        alpha_0[i] = 0.1;
                        alpha_1[i] = 1.0;
                }
                priors.push_back(alpha_0);
                priors.push_back(alpha_1);

                transition.push_back(phylotree_hmm_t<code_t, alphabet_size>::vector_t(2, 0.0));
                transition.push_back(phylotree_hmm_t<code_t, alphabet_size>::vector_t(2, 0.0));

                transition[0][0] = 0.95; // 0 -> 0
                transition[0][1] = 0.05; // 0 -> 1

                transition[1][0] = 0.05; // 1 -> 0
                transition[1][1] = 0.95; // 1 -> 1
        }
} options_t;

static options_t options;

// Main
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] TREE FASTA_ALIGNMENT\n\n", pname);
        (void)fprintf(fp,
                      "\n"
                      "Options:\n"
                      "             -a VECTOR       - pseudo count\n"
                      "             -d DIMENSION    - dimension\n"
                      "             -t MATRIX       - transition matrix\n"
                      "\n"
                      "             -v              - be verbose\n"
                      "   --help                    - print help and exit\n"
                      "   --version                 - print version information and exit\n\n");
}

static
void wrong_usage(const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `phylotree-hmm-test --help' for more information.\n");

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

static
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

static
void run_hmm(const char* file_tree, const char* file_alignment)
{
        size_t dimension = options.dimension;

        /* phylogenetic tree */
        pt_root_t* pt_root = parse_tree_file(file_tree);
        pt_root->init(alphabet_size);

        /* alignment */
        alignment_t<code_t> alignment(file_alignment, pt_root);

        /* uniform distribution on the initial state */
        phylotree_hmm_t<code_t, alphabet_size>::vector_t px_0(dimension, 1.0/(double)dimension);

        phylotree_hmm_t<code_t, alphabet_size> hmm(pt_root, px_0, options.transition, options.priors);
        hmm.run(alignment);

        for (size_t i = 0; i < hmm.size(); i++) {
                cout << setprecision(8)
                     << fixed
                     << hmm[i][0] << " "
                     << hmm[i][1]
                     << endl;
        }
}

void init_options(const string& alpha, const string& transition)
{
        vector<string> tmp;
        vector<string> tmp_;

        if (alpha != "") {
                options.priors = phylotree_hmm_t<code_t, alphabet_size>::priors_t();
                tmp = token(strip(alpha), ' ');
                for (size_t i = 0; i < tmp.size(); i++) {
                        if (tmp[i] == "") {
                                continue;
                        }
                        string str = tmp[i];
                        exponent_t<code_t, alphabet_size> alpha;
                        for (size_t i = 0; i < alphabet_size; i++) {
                                alpha[i] = atof(str.c_str());
                        }
                        options.priors.push_back(alpha);
                }
                if (options.priors.size() != options.dimension) {
                        cerr << "Dimensions do not match."
                             << endl;
                        exit(EXIT_FAILURE);
                }
        }
        if (options.verbose) {
                for (size_t i = 0; i < options.dimension; i++) {
                        cerr << "alpha[" << i << "]: ";
                        for (size_t j = 0; j < alphabet_size; j++) {
                                cerr << options.priors[i][j]
                                     << " ";
                        }
                        cerr << endl;
                }
        }
        if (transition != "") {
                options.transition = phylotree_hmm_t<code_t, alphabet_size>::matrix_t();
                for (size_t i = 0; i < options.dimension; i++) {
                        options.transition.push_back(
                                phylotree_hmm_t<code_t, alphabet_size>::vector_t(options.dimension, 0.0));
                }
                tmp = token(strip(transition), ';');
                for (size_t i = 0; i < tmp.size(); i++) {
                        tmp_ = token(tmp[i], ' ');
                        for (size_t j = 0; j < tmp_.size(); j++) {
                                string str = tmp_[j];
                                if (i >= options.dimension || j >= options.dimension) {
                                        cerr << "Dimensions do not match."
                                             << endl;
                                        exit(EXIT_FAILURE);
                                }
                                options.transition[i][j] = atof(str.c_str());
                        }
                }
        }
        if (options.verbose) {
                cerr << "transition matrix:"
                     << endl;
                for (size_t i = 0; i < options.transition.size(); i++) {
                        for (size_t j = 0; j < options.transition[0].size(); j++) {
                                cerr << setprecision(3)
                                     << fixed
                                     << options.transition[i][j]
                                     << " ";
                        }
                        cerr << endl;
                }
        }
}

int main(int argc, char *argv[])
{
        const char* file_tree;
        const char* file_alignment;

        string alpha;
        string transition;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 's' }
                };

                c = getopt_long(argc, argv, "a:d:t:vh",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        alpha = string(optarg);
                        break;
                case 'd':
                        options.dimension = atoi(optarg);
                        break;
                case 't':
                        transition = string(optarg);
                        break;
                case 'v':
                        options.verbose = true;
                        break;
                case 'h':
                        print_usage(argv[0], stdout);
                        exit(EXIT_SUCCESS);
                case 's':
                        print_version(stdout);
                        exit(EXIT_SUCCESS);
                default:
                        wrong_usage(NULL);
                        exit(EXIT_FAILURE);
                }
        }
        if(optind+2 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }
        init_options(alpha, transition);

        file_tree      = argv[optind+0];
        file_alignment = argv[optind+1];

        run_hmm(file_tree, file_alignment);

        return 0;
}
