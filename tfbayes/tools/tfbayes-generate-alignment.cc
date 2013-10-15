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
#include <tfbayes/uipac/alphabet.hh>
#include <tfbayes/utility/strtools.hh>
#include <tfbayes/exception/exception.h>

#define alphabet_size 5

using namespace std;
using boost::format;

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        vector<double> alpha; // pseudocounts for the motifs
        vector<double> beta;  // pseudocounts for the background
        double d;
        double lambda;
        size_t n;
        size_t m;
        string format;
        _options_t()
                : alpha(alphabet_size, 0.2),
                  beta (alphabet_size, 2.0),
                  d(1.0),
                  lambda(0.01),
                  n(100),
                  m(10),
                  format("pretty")
                { }
} options_t;

static options_t options;

// Usage
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] MODEL TREE\n\n", pname);
        (void)fprintf(fp,
                      "Models:\n"
                      "             simple           - all columns are generated with the same pseudocounts"
                      "             tfbs             - enrich alignment with certain patterns"
                      "\n"
                      "Options:\n"
                      "             -a F:F:F:F:F     - pseudo counts (five floats separated by a colon)\n"
                      "             -b F:F:F:F:F     - background pseudo counts\n"
                      "             -d CONCENTRATION - dirichlet process concentration parameter\n"
                      "             -l WEIGHT        - weight of the foreground\n"
                      "             -n LENGTH        - length of the alignment\n"
                      "             -m LENGTH        - tfbs/pattern length\n"
                      "             --format=FORMAT  - output format (pretty, fasta)\n"
                      "\n"
                      "   --help                     - print help and exit\n"
                      "   --version                  - print version information and exit\n\n");
}

static
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

// Alignments
////////////////////////////////////////////////////////////////////////////////

void print_alignment(const alignment_t<>& alignment)
{
        if (options.format == "pretty") {
                cout << print_alignment_pretty(alignment) << endl;
        }
        else if (options.format == "fasta") {
                cout << print_alignment_fasta(alignment) << endl;
        }
        else {
                wrong_usage("Unknown output format.");
                exit(EXIT_FAILURE);
        }
}

class pattern_t : public matrix<double> {
public:
        pattern_t()
                : matrix<double>(options.m, alphabet_size, 0.0)
                { }
        pattern_t(const vector<double>& pseudocounts, gsl_rng * r)
                : matrix<double>(options.m, 0) {
                // generate new pattern
                for (size_t i = 0; i < options.m; i++) {
                        operator[](i) = dirichlet_sample<alphabet_size>(pseudocounts, r);
                }
        }
};

class dirichlet_process_t {
public:
        dirichlet_process_t(double alpha, const vector<double>& pseudocounts, gsl_rng * r)
                : alpha(alpha), n(0), pseudocounts(pseudocounts), r(r) {
        }
        const pattern_t& operator()() {
                size_t c   = occurrences.size();
                double sum = 0.0;
                double ra  = (double)rand()/RAND_MAX;

                // draw from an existing cluster
                for (size_t i = 0; i < c; i++) {
                        sum += occurrences[i]/(n + alpha);
                        if (ra <= sum) {
                                occurrences[i] += 1.0;
                                n += 1.0;
                                return patterns[i];
                        }
                }
                // generate new pattern
                patterns   .push_back(pattern_t(pseudocounts, r));
                occurrences.push_back(1);
                n += 1.0;

                return patterns[c];
        }
        friend ostream& operator<<(ostream& o, const dirichlet_process_t& dirichlet_process);

protected:
        double alpha;
        double n;
        vector<double> pseudocounts;
        vector<double> occurrences;
        vector<pattern_t> patterns;
        gsl_rng * r;
};

ostream& operator<<(ostream& o, const pattern_t& pattern)
{
        for (size_t j = 0; j < pattern[0].size(); j++) {
                o << format("%c: ") % nucleotide_alphabet_t().decode(j);
                for (size_t i = 0; i < pattern.size(); i++) {
                        o << format("%5f ") % pattern[i][j];
                }
                o << endl;
        }
        return o;
}

ostream& operator<<(ostream& o, const dirichlet_process_t& dirichlet_process)
{
        o << "Cluster: " << dirichlet_process.occurrences.size() << endl
          << "Total N: " << dirichlet_process.n                  << endl
          << endl;

        for (size_t k = 0; k < dirichlet_process.occurrences.size(); k++) {
                const pattern_t& pattern = dirichlet_process.patterns[k];
                o << format("Cluster %d (occurred %d times)\n") % k % dirichlet_process.occurrences[k];
                o << "-----------------------------------------------------" << endl;
                o << pattern << endl;
        }
        return o;
}

void insert_observations(alignment_t<>& alignment, size_t i, const vector<alphabet_code_t>& observations)
{
        for (size_t k = 0; k < observations.size(); k++) {
                alignment[alignment_index_t(k, i)] = observations[k];
        }
}

void generate_tfbs_alignment(const pt_root_t& pt_root, gsl_rng * r)
{
        dirichlet_process_t dirichlet_process(options.d, options.alpha, r);
        vector<double> stationary;
        vector<alphabet_code_t> observations;
        // alignment
        alignment_t<> alignment(options.n, pt_root);

        // generate
        for (size_t i = 0; i < options.n; i++) {
                // generate foreground
                if ((double)rand()/RAND_MAX <= options.lambda) {
                        const pattern_t& pattern = dirichlet_process();
                        for (size_t j = 0; j < pattern.size() && i+j < options.n; j++) {
                                observations = pt_generate_observations<alphabet_size, alphabet_code_t>(pt_root, pattern[j]);
                                insert_observations(alignment, i+j, observations);
                        }
                        i += pattern.size() - 1;
                }
                // generate background
                else {
                        stationary   = dirichlet_sample<alphabet_size>(options.beta, r);
                        observations = pt_generate_observations<alphabet_size, alphabet_code_t>(pt_root, stationary);
                        insert_observations(alignment, i, observations);
                }
        }
        // print cluster
        cerr << dirichlet_process << endl;
        // print result
        print_alignment(alignment);
}

void generate_simple_alignment(const pt_root_t& pt_root, gsl_rng * r)
{
        // alignment
        alignment_t<> alignment(options.n, pt_root);

        // generate
        for (size_t i = 0; i < options.n; i++) {
                vector<double         > stationary   =
                        dirichlet_sample<alphabet_size>(options.alpha, r);
                vector<alphabet_code_t> observations =
                        pt_generate_observations<alphabet_size, alphabet_code_t>(pt_root, stationary);

                alignment[i] = observations;
        }
        // print result
        print_alignment(alignment);
}

// Main
////////////////////////////////////////////////////////////////////////////////

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

pt_root_t parse_tree_file(const char* file_tree)
{
        FILE* yyin = fopen(file_tree, "r");
        if (yyin == NULL) {
                std_err(PERR, "Could not open phylogenetic tree");
        }
        list<pt_root_t> tree_list = parse_tree_list(yyin);
        assert(tree_list.size() == 1);
        pt_root_t pt_root = tree_list.front();
        fclose(yyin);

        return pt_root;
}

void generate_alignment(const string& model, const char* treefile)
{
        init();

        // gsl random number generator
        gsl_rng_env_setup();
        const gsl_rng_type * T = gsl_rng_default;
        gsl_rng * r = gsl_rng_alloc(T);

        // parse tree
        pt_root_t pt_root = parse_tree_file(treefile);

        // switch command
        if (model == "simple") {
                generate_simple_alignment(pt_root, r);
        }
        else if (model == "tfbs") {
                generate_tfbs_alignment(pt_root, r);
        }
        else {
                wrong_usage("Invalid statistical model.");
                exit(EXIT_FAILURE);
        }

        // clean up
        gsl_rng_free (r);
}

int main(int argc, char *argv[])
{
        vector<string> tokens;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "format",          1, 0, 'f' },
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'v' }
                };

                c = getopt_long(argc, argv, "a:b:d:l:n:hv",
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
                case 'b':
                        tokens = token(string(optarg), ':');
                        if (tokens.size() != alphabet_size) {
                                wrong_usage(NULL);
                                exit(EXIT_FAILURE);
                        }
                        for (size_t i = 0; i < alphabet_size; i++) {
                                options.beta[i] = atof(tokens[i].c_str());
                        }
                        break;
                case 'f':
                        options.format = string(optarg);
                        break;
                case 'd':
                        options.d = atof(optarg);
                        break;
                case 'l':
                        options.lambda = atof(optarg);
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
        if(optind+2 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }
        generate_alignment(string(argv[optind]), argv[optind+1]);

        return 0;
}
