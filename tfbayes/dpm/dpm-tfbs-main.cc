/* Copyright (C) 2011 Philipp Benner
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

#include <ctime>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include <getopt.h>

#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-sampler.hh>
#include <tfbayes/dpm/dpm-sampling-history.hh>
#include <tfbayes/fasta/fasta.hh>
#include <tfbayes/utility/strtools.hh>

using namespace std;

typedef struct _options_t {
        size_t samples;
        size_t burnin;
        size_t foreground_length_min;
        size_t foreground_length_max;
        double alpha;
        double discount;
        double lambda;
        const char* process_prior;
        const char* background_model;
        double background_alpha;
        size_t background_context;
        size_t population_size;
        string save;
        _options_t()
                : samples(1000),
                  burnin(100),
                  foreground_length_min( 8),
                  foreground_length_max(10),
                  alpha(0.05),
                  discount(0.0),
                  lambda(0.01),
                  process_prior("pitman-yor process"),
                  background_model("dirichlet-mixture"),
                  background_alpha(1),
                  background_context(2),
                  population_size(1),
                  save()
                { }
} options_t;

ostream&
operator<<(std::ostream& o, const _options_t& options) {
        o << "Options:"              << endl
          << "-> samples             = " << options.samples             << endl
          << "-> burnin              = " << options.burnin              << endl
          << "-> alpha               = " << options.alpha               << endl
          << "-> discount            = " << options.discount            << endl
          << "-> lambda              = " << options.lambda              << endl
          << "-> foreground length min     = " << options.foreground_length_min     << endl
          << "-> foreground length max     = " << options.foreground_length_max     << endl
          << "-> process prior       = " << options.process_prior       << endl
          << "-> background model    = " << options.background_model    << endl
          << "-> background_alpha    = " << options.background_alpha    << endl
          << "-> background_context  = " << options.background_context  << endl
          << "-> population_size     = " << options.population_size     << endl
          << "-> save                = " << options.save                << endl;
        return o;
}

static options_t options;

static
void print_usage(char *pname, FILE *fp)
{
	(void)fprintf(fp,
                      "\nUsage: %s [OPTION]... PHYLOGENETIC_DATA FASTA_ALIGNMENT\n\n", pname);
	(void)fprintf(fp,
                      "Options:\n"
                      "   --alpha=ALPHA             - alpha parameter for the dirichlet process\n"
                      "   --lambda=LAMBDA           - lambda mixture weight\n"
                      "   --tfbs-length=TFBS_LENGTH - length of the tfbs\n"
                      "   --population-size=N       - number of parallel samplers\n"
                      "\n"
                      "   --samples=SAMPLES:BURN_IN - number of samples\n"
                      "   --save=FILE_NAME          - save posterior to file\n"
                      "\n"
                      "   --help	            - print help and exit\n"
                      "   --version	            - print version information and exit\n\n");
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

static
void wrong_usage(const char *msg)
{

	if(msg != NULL) {
		(void)fprintf(stderr, "%s\n", msg);
	}
	(void)fprintf(stderr,
		      "Try `dpm-tfbs-debug --help' for more information.\n");

	exit(EXIT_FAILURE);

}

ostream& operator<< (ostream& o, const product_dirichlet_t& pd) {
        for (size_t j = 0; j < pd.alpha[0].size() - 1; j++) {
                o << "\t";
                for (size_t i = 0; i < pd.alpha.size(); i++) {
                        o << pd.alpha[i][j] + pd.counts[i][j] << " ";
                }
                o << endl;
        }

        return o;
}

static
void run_dpm(const char* phylogenetic_data_file, const char* fasta_alignment_file)
{
        tfbs_options_t tfbs_options;

        // background alpha
        tfbs_options.background_alpha.push_back(vector<double>(data_tfbs_t::alphabet_size, options.background_alpha));

        // baseline prior
        tfbs_options.baseline_priors.push_back(matrix<double>());
        tfbs_options.baseline_priors.begin()->push_back(vector<double>(data_tfbs_t::alphabet_size, 0.2));

        // tfbs options
        tfbs_options.phylogenetic_file   = phylogenetic_data_file;
        tfbs_options.alignment_file      = fasta_alignment_file;
        tfbs_options.alpha               = options.alpha;
        tfbs_options.lambda              = options.lambda;
        tfbs_options.discount            = options.discount;
        tfbs_options.population_size     = options.population_size;
        tfbs_options.process_prior       = options.process_prior;
        tfbs_options.background_model    = options.background_model;
        tfbs_options.background_gamma    = vector<double>(2,1);
        tfbs_options.background_context  = options.background_context;
        tfbs_options.baseline_weights    = vector<double>(1,1);
        tfbs_options.baseline_names.push_back("baseline_default");
        tfbs_options.block_samples       = false;
        tfbs_options.block_samples_period= 1;
        tfbs_options.metropolis_proposals= 4;
        tfbs_options.optimize            = false;
        tfbs_options.optimize_period     = 1;
        tfbs_options.initial_temperature = 1.0;
        tfbs_options.threads             = 1;
        tfbs_options.verbose             = true;
        tfbs_options.baseline_lengths.push_back(vector<double>());
        for (size_t i = options.foreground_length_min; i <= options.foreground_length_max; i++) {
                tfbs_options.baseline_lengths[0].push_back(i);
        }

        // create data, dpm, and sampler objects
        dpm_tfbs_pmcmc_t pmcmc(tfbs_options);

        // execute the sampler
        pmcmc(options.samples, options.burnin);

        // save result
        pmcmc.save(options.save);
}

int main(int argc, char *argv[])
{
        char *phylogenetic_data_file;
	char *fasta_alignment_file;
        vector<string> tokens;

	if(argc == 1) {
		wrong_usage("Too few arguments.");
		exit(EXIT_FAILURE);
	}

	for(;;) {
		int c, option_index = 0;
		static struct option long_options[] = {
                        { "alpha",           1, 0, 'a' },
                        { "lambda",          1, 0, 'l' },
                        { "context",         1, 0, 'c' },
                        { "samples",         1, 0, 's' },
                        { "tfbs-length",     1, 0, 't' },
                        { "population-size", 1, 0, 'p' },
                        { "save",            1, 0, 'e' },
			{ "help",	     0, 0, 'h' },
			{ "version",	     0, 0, 'v' }
		};

		c = getopt_long(argc, argv, "d:",
				long_options, &option_index);

		if(c == -1) {
			break;
		}

		switch(c) {
                case 'a':
                        options.alpha = atof(optarg);
                        break;
                case 'd':
                        options.discount = atof(optarg);
                        if (options.discount < 0 || options.discount >= 1) {
                                wrong_usage(NULL);
                        }
                        break;
                case 'l':
                        options.lambda = atof(optarg);
                        break;
                case 'c':
                        options.background_context = atoi(optarg);
                        break;
                case 'e':
                        options.save = string(optarg);
                        break;
                case 's':
                        if (token(optarg, ':').size() != 2) {
                                wrong_usage(NULL);
                        }
                        options.samples = atoi(token(optarg, ':')[0].c_str());
                        options.burnin  = atoi(token(optarg, ':')[1].c_str());
                        break;
                case 't':
                        tokens = token(string(optarg), ':');
                        if (tokens.size() == 1) {
                                options.foreground_length_min = atoi(tokens[0].c_str());
                                options.foreground_length_max = atoi(tokens[0].c_str());
                        }
                        else if (tokens.size() == 2) {
                                options.foreground_length_min = atoi(tokens[0].c_str());
                                options.foreground_length_max = atoi(tokens[1].c_str());
                        }
                        else {
                                wrong_usage(NULL);
                        }
                        break;
                case 'p':
                        options.population_size = atoi(optarg);
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

        cout << options << endl;

	phylogenetic_data_file = argv[optind+0];
        fasta_alignment_file   = argv[optind+1];

        run_dpm(phylogenetic_data_file, fasta_alignment_file);

        return 0;
}
