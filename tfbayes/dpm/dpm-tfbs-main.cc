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
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <ctime>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include <getopt.h>

#include <init.hh>
#include <dpm-tfbs.hh>
#include <pmcmc.hh>

#include <tfbayes/exception.h>
#include <tfbayes/fasta.hh>

using namespace std;

typedef struct _options_t {
        size_t samples;
        size_t burnin;
        size_t tfbs_length;
        double alpha;
        double d;
        double lambda;
        size_t population_size;
        string save;
        _options_t()
                : samples(1000),
                  burnin(100),
                  tfbs_length(10),
                  alpha(0.05),
                  d(0.0),
                  lambda(0.01),
                  population_size(1),
                  save()
                { }
} options_t;

ostream&
operator<<(std::ostream& o, const _options_t& options) {
        o << "Options:"              << endl
          << "-> samples         = " << options.samples         << endl
          << "-> burnin          = " << options.burnin          << endl
          << "-> tfbs_length     = " << options.tfbs_length     << endl
          << "-> alpha           = " << options.alpha           << endl
          << "-> d               = " << options.d               << endl
          << "-> lambda          = " << options.lambda          << endl
          << "-> population_size = " << options.population_size << endl
          << "-> save            = " << options.save            << endl;
        return o;
}

static options_t options;

static
void print_usage(char *pname, FILE *fp)
{
	(void)fprintf(fp,
                      "\nUsage: %s [OPTION]... FILE\n\n", pname);
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
		      "Try `dpm-tfbs --help' for more information.\n");

	exit(EXIT_FAILURE);

}

static
char * readfile(const char* file_name, vector<string>& sequences)
{
        FastaParser parser(file_name);

        ifstream file(file_name);
        string line;

        size_t i = 0;
        while ((line = parser.read_sequence()) != "") {
                size_t read = line.size();
                sequences.push_back("");
                size_t pos = 0;
                for (size_t j = 0; j < (size_t)read && line[j] != '\n'; j++) {
                        if (is_nucleotide_or_masked(line[j])) {
                                sequences[i].append(1,line[j]);
                                pos++;
                        }
                }
                i++;
        }

        return NULL;
}

ostream& operator<< (ostream& o, const ProductDirichlet& pd) {
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
void save_motifs(ostream& file, const DpmTfbs& dpm)
{
        const ClusterManager& cm = dpm.clustermanager();

        for (ClusterManager::const_iterator it = cm.begin();
             it != cm.end(); it++) {
                if ((*it)->cluster_tag() == 0) {
                        file << "cluster_bg" << " =" << endl;
                        file << static_cast<const ProductDirichlet&>((*it)->model());
                }
        }
}

static
void save_result(ostream& file, Sampler& sampler)
{
        const posterior_t& posterior      = sampler.posterior();
        const sampling_history_t& history = sampler.sampling_history();
        file.setf(ios::showpoint);

        file << "[Result]" << endl;
        file << "posterior =" << endl;
        for (size_t i = 0; i < posterior.probabilities.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < posterior.probabilities[i].size(); j++) {
                        file << (float)posterior.probabilities[i][j] << " ";
                }
                file << endl;
        }
        file << "components =" << endl;
        for (size_t i = 0; i < history.components.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.components[i].size(); j++) {
                        file << history.components[i][j] << " ";
                }
                file << endl;
        }
        file << "switches =" << endl;
        for (size_t i = 0; i < history.switches.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.switches[i].size(); j++) {
                        file << history.switches[i][j] << " ";
                }
                file << endl;
        }
        file << "likelihood =" << endl;
        for (size_t i = 0; i < history.likelihood.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < history.likelihood[i].size(); j++) {
                        file << history.likelihood[i][j] << " ";
                }
                file << endl;
        }
        file << "graph = ";
        for (Graph::const_iterator it = posterior.graph.begin();
             it != posterior.graph.end(); it++) {
                file << (*it).first.index1 << "-"
                     << (*it).first.index2 << "="
                     << static_cast<double>((*it).second)/static_cast<double>(sampler.sampling_steps()) << " ";
        }
        file << endl;
}

static
void run_dpm(const char* file_name)
{
        vector<string> sequences;
        vector<string> sequences_comp;

        // read sequences
        readfile(file_name, sequences);
        sequences_comp = complement(sequences);

        // baseline
        vector<double> baseline_weights(1,1);
        gsl_matrix *baseline_priors[2];
        baseline_priors[0] = gsl_matrix_alloc(options.tfbs_length, 4);
        baseline_priors[1] = NULL;
        for (size_t i = 0; i < options.tfbs_length; i++) {
                for (size_t j = 0; j < 4; j++) {
                        gsl_matrix_set(baseline_priors[0], i, j, 1.0);
                }
        }

        // create data, dpm, and sampler objects
        data_tfbs_t& data = *new data_tfbs_t(sequences, options.tfbs_length);
        data_tfbs_t& data_comp = *new data_tfbs_t(sequences_comp, options.tfbs_length);
        DpmTfbs& gdpm = *new DpmTfbs(options.alpha, options.d, options.lambda, options.tfbs_length, data, data_comp, baseline_weights, baseline_priors);
        GibbsSampler& sampler = *new GibbsSampler(gdpm, data);
        PopulationMCMC& pmcmc = *new PopulationMCMC(sampler, options.population_size);

        // execute the sampler
        pmcmc.sample(options.samples, options.burnin);

        // save result
        if (options.save == "") {
                save_result(cout, pmcmc);
                save_motifs(cout, gdpm);
        }
        else {
                ofstream file;
                file.open(options.save.c_str());
                save_result(file, pmcmc);
                save_motifs(file, gdpm);
                file.close();
        }

        // free memory
        free(&data);
        free(&pmcmc);
}

static
vector<string> token(const string& str, char t) {
        string token;
        vector<string> tokens;
        istringstream iss(str);
        while (getline(iss, token, t)) {
                tokens.push_back(token);
        }
        return tokens;
}

int main(int argc, char *argv[])
{
        __dpm_init__();

	char *file_name;

	if(argc == 1) {
		wrong_usage("Too few arguments.");
		exit(EXIT_FAILURE);
	}

	for(;;) {
		int c, option_index = 0;
		static struct option long_options[] = {
                        { "alpha",           1, 0, 'a' },
                        { "lambda",          1, 0, 'l' },
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
                        options.d = atof(optarg);
                        if (options.d < 0 || options.d >= 1) {
                                wrong_usage(NULL);
                        }
                        break;
                case 'l':
                        options.lambda = atof(optarg);
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
                        options.tfbs_length = atoi(optarg);
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
        cout << options << endl;

	file_name = argv[optind];

        run_dpm(file_name);

        return 0;
}
