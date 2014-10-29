/* Copyright (C) 2014 Philipp Benner
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

#include <algorithm>
#include <sstream>
#include <getopt.h>
#include <string.h>
#include <math.h>

#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/dpm/component-model.hh>
#include <tfbayes/dpm/dpm-tfbs-options.hh>
#include <tfbayes/uipac/alphabet.hh>

using namespace std;

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        double alpha;
        size_t context;
        bool   verbose;
        _options_t()
                : alpha(1.0),
                  context(1),
                  verbose(false)
                { }
} options_t;

static options_t options;

// Usage
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] FASTA\n\n", pname);
        (void)fprintf(fp,
                      "\n"
                      "Options:\n"
                      "             -a double       - pseudocounts\n"
                      "             -n integer      - maximum context\n"
                      "\n"
                      "             -v              - be verbose\n"
                      "   --help                    - print help and exit\n"
                      "   --version                 - print version information and exit\n\n");
}

void wrong_usage(const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `mcm-test --help' for more information.\n");

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

// Test functions
////////////////////////////////////////////////////////////////////////////////

static sequence_data_t<data_tfbs_t::code_t>
read_data(const string& fasta_file, alphabet_t& alphabet) {
        FastaParser parser(fasta_file);
        sequence_data_t<data_tfbs_t::code_t> data;

        while (parser) {
                string line = parser();
                if (options.verbose) {
                        if (parser.description().size() > 0) {
                                cerr << "Reading sequence `"
                                     << parser.description()[0]
                                     << "'." << endl;
                        }
                        else {
                                cerr << "Reading sequence without description."
                                     << endl;
                        }
                }
                sequence_t<alphabet_code_t> sequence(line, alphabet);
                vector<data_tfbs_t::code_t> tmp(sequence.size());
                // loop over the sequence
                for (size_t i = 0; i < sequence.size(); i++) {
                        tmp[i][sequence[i]] = 1;
                }
                data.push_back(tmp);
        }

        return data;
}

static void
_main_(const string& fasta_file) {
        alphabet_t alphabet = nucleotide_alphabet_t();
        sequence_data_t<data_tfbs_t::code_t> data
                = read_data(fasta_file, alphabet);
        sequence_data_t<cluster_tag_t> cluster_assignments(data.sizes(), 0);

        // the alignment should contain at least two sequences, one
        // for learning the model, and another one for predictions
        assert(data.size() > 1);

        tfbs_options_t tfbs_options;
        tfbs_options.background_alpha   = matrix<double>(1,1,options.alpha);
        tfbs_options.background_context = options.context;
        tfbs_options.background_weights = "entropy";

        markov_chain_mixture_t model(alphabet.size(), tfbs_options, data,
                                     cluster_assignments, 1);

        for (size_t i = 0; i < data.size()-1; i++) {
                range_t range(seq_index_t(i, 0), data[i].size());
        }
        range_t range(seq_index_t(data.size()-1, 0), data[data.size()-1].size());
        cout << "Log predictive: "
             << model.log_predictive(range)
             << endl;
}

int main(int argc, char *argv[])
{
        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",              0, 0, 'h' },
                        { "version",           0, 0, 'q' }
                };

                c = getopt_long(argc, argv, "a:n:v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        options.alpha = atof(optarg);
                        break;
                case 'n':
                        options.context = atoi(optarg);
                        break;
                case 'v':
                        options.verbose = true;
                        break;
                case 'h':
                        print_usage(argv[0], stdout);
                        exit(EXIT_SUCCESS);
                case 'q':
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
        _main_(argv[optind]);

        return 0;
}
