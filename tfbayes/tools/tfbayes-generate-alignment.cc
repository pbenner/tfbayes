/* Copyright (C) 2013 Philipp Benner
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
#include <iomanip>
#include <vector>
#include <sys/time.h>
#include <getopt.h>
#include <cstdlib>

#include <boost/format.hpp>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-generate-observations.hh>
#include <tfbayes/uipac/code.hh>
#include <tfbayes/utility/strtools.hh>
#include <tfbayes/exception/exception.h>

#define alphabet_size 5
typedef short code_t;

using namespace std;
using boost::format;

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

void print_alignment(const pt_root_t* root, const alignment_t<code_t>& alignment)
{
        for (size_t i = 0; i < alignment.n_species(); i++) {
                cout << format("%10s: ")
                        % root->leaves[i]->name;

                for (size_t j = 0; j < alignment.length(); j++) {
                        cout << decode_nucleotide(alignment[i][j]);
                }
                cout << endl;
        }
}

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        vector<double> alpha;
        size_t n;
        _options_t()
                : alpha(alphabet_size, 0.2),
                  n(100)
                { }
} options_t;

static options_t options;

// Main
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] TREE\n\n", pname);
        (void)fprintf(fp,
                      "Methods: gradient-ascent, metropolis-hastings\n"
                      "\n"
                      "Options:\n"
                      "             -a F:F:F:F:F    - pseudo counts (five floats separated by a colon)\n"
                      "             -n LENGTH       - length of the alignment\n"
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
                      "Try `tfbayes-generate-alignment --help' for more information.\n");

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

pt_root_t* parse_tree_file(const char* file_tree)
{
        FILE* yyin = fopen(file_tree, "r");
        if (yyin == NULL) {
                std_err(PERR, "Could not open phylogenetic tree");
        }
        list<pt_root_t*> tree_list = parse_tree_list(yyin);
        assert(tree_list.size() == 1);
        pt_root_t* pt_root = tree_list.front();
        fclose(yyin);

        return pt_root;
}

void generate_alignment(const char* treefile)
{
        init();

        // gsl random number generator
        gsl_rng_env_setup();
        const gsl_rng_type * T = gsl_rng_default;
        gsl_rng * r = gsl_rng_alloc(T);

        // parse tree
        pt_root_t* pt_root = parse_tree_file(treefile);

        // alignment
        alignment_t<code_t> alignment(options.n, -1, pt_root);

        // generate
        for (size_t i = 0; i < options.n; i++) {
                vector<double> stationary   = dirichlet_sample<alphabet_size>(options.alpha, r);
                vector<code_t> observations = pt_generate_observations<code_t, alphabet_size>(pt_root, stationary);

                for (size_t k = 0; k < observations.size(); k++) {
                        alignment[k][i] = observations[k];
                }
        }
        // print result
        print_alignment(pt_root, alignment);

        // clean up
        pt_root->destroy();
        gsl_rng_free (r);
}

int main(int argc, char *argv[])
{
        vector<string> tokens;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'v' }
                };

                c = getopt_long(argc, argv, "a:b:e:m:n:r:l:p:s:j:hv",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        tokens = token(string(optarg), ':');
                        if (tokens.size() != alphabet_size) {
                                wrong_usage(NULL);
                                exit(EXIT_FAILURE);
                        }
                        for (size_t i = 0; i < alphabet_size; i++) {
                                options.alpha[i] = atof(tokens[i].c_str());
                        }
                        break;
                case 'n':
                        options.n = atoi(optarg);
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
        if(optind+1 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }
        generate_alignment(argv[optind]);

        return 0;
}
