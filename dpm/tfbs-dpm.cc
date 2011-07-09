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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include <sampler.hh>
#include <init.hh>
#include <tfbayes/exception.h>

using namespace std;

typedef struct _options_t {
        size_t samples;
        size_t burnin;
        size_t tfbs_length;
        double alpha;
        double lambda;
        string save;
        _options_t()
                : samples(1000),
                  burnin(100),
                  tfbs_length(10),
                  alpha(0.05),
                  lambda(0.01),
                  save()
                { }
} options_t;

ostream&
operator<<(std::ostream& o, const _options_t& options) {
        o << "Options:"          << endl
          << "-> samples     = " << options.samples     << endl
          << "-> burnin      = " << options.burnin      << endl
          << "-> tfbs_length = " << options.tfbs_length << endl
          << "-> alpha       = " << options.alpha       << endl
          << "-> lambda      = " << options.lambda      << endl
          << "-> save        = " << options.save        << endl;
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
                      "\n"
                      "   --samples=SAMPLES:BURN_IN - number of samples\n"
                      "   --save=FILE_NAME          - save posterior to file\n"
                      "\n"
                      "   --help	            - print help and exit\n"
                      "   --version	            - print version information and exit\n\n");

	return;

}

static
void print_version(FILE *fp)
{

	(void)fprintf(fp,
                      "This is free software, and you are welcome to redistribute it\n"
                      "under certain conditions; see the source for copying conditions.\n"
                      "There is NO warranty; not even for MERCHANTABILITY or FITNESS\n"
                      "FOR A PARTICULAR PURPOSE.\n\n");

	return;

}

static
void wrong_usage(const char *msg)
{

	if(msg != NULL) {
		(void)fprintf(stderr, "%s\n", msg);
	}
	(void)fprintf(stderr,
		      "Try `tfbs-dpm --help' for more information.\n");

	exit(EXIT_FAILURE);

}

static inline
int is_nucleotide(char S) {
        if (S == 'A' || S == 'C' || S == 'G' || S == 'T' ||
            S == 'a' || S == 'c' || S == 'g' || S == 't') {
                return 1;
        }
        else {
                return 0;
        }
}

static
void get_max_length(const char* file_name, size_t *lines, size_t *max_len)
{
        string line;
        ifstream file(file_name);

        *max_len = 0;
        *lines   = 0;

        while (getline(file, line)) {
                size_t read = line.size();
                if (line[read-1] == '\n') {
                        read--;
                }
                if ((size_t)read > *max_len) {
                        *max_len = (size_t)read;
                }
                if (read > 0) {
                        (*lines)++;
                }
        }

}

static
char * readfile(const char* file_name, char *sequences[])
{
        ifstream file(file_name);
        string line;

        size_t i = 0;
        while (getline(file, line)) {
                size_t read = line.size();
                if (read > 1) {
                        size_t pos = 0;
                        for (size_t j = 0; j < (size_t)read && line[j] != '\n'; j++) {
                                if (is_nucleotide(line[j])) {
                                        sequences[i][pos] = line[j];
                                        pos++;
                                }
                        }
                        i++;
                }
        }

        return NULL;
}

static
char ** alloc_sequences(size_t n, size_t m) {
        char **sequences = (char **)malloc((n+1)*sizeof(char *));
        size_t i;

        for (i = 0; i < n; i++) {
                sequences[i] = (char *)calloc(m, sizeof(char));
        }
        sequences[i] = NULL;

        return sequences;
}

static
void free_sequences(char **sequences)
{
        for (size_t i = 0; sequences[i] != NULL; i++) {
                free(sequences[i]);
        }
}

static
void run_dpm(const char* file_name)
{
        size_t lines, max_len;
        char **sequences;

        get_max_length(file_name, &lines, &max_len);
        sequences = alloc_sequences(lines, max_len);

        readfile(file_name, sequences);

        Data* data = new Data(lines, sequences);
        DPM*  gdpm = new DPM(options.alpha, options.lambda, *data);
        GibbsSampler* sampler = new GibbsSampler(*gdpm, *data);

        sampler->sample(options.samples, options.burnin);
        const vector<vector<double> >& posterior = gdpm->posterior();

        if (options.save == "") {
                for (size_t i = 0; i < posterior.size(); i++) {
                        for (size_t j = 0; j < posterior[i].size(); j++) {
                                printf("%f ", (float)posterior[i][j]);
                        }
                        cout << endl;
                }
        }
        else {
                ofstream file;
                file.open(options.save.c_str());
                file.setf(ios::showpoint); 
                for (size_t i = 0; i < posterior.size(); i++) {
                        for (size_t j = 0; j < posterior[i].size(); j++) {
                                file << (float)posterior[i][j] << " ";
                        }
                        file << endl;
                }
                file.close();
        }

        free(data);
        free(gdpm);
        free(sampler);
        free_sequences(sequences);
}

static
string token(const string& str, size_t i) {
        string token;
        vector<string> tokens;
        istringstream iss(str);
        while (getline(iss, token, ':')) {
                tokens.push_back(token);
        }
        return tokens[i];
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
                        { "alpha",       1, 0, 'a' },
                        { "lambda",      1, 0, 'l' },
                        { "samples",     1, 0, 's' },
                        { "tfbs-length", 1, 0, 't' },
                        { "save",        1, 0, 'e' },
			{ "help",	 0, 0, 'h' },
			{ "version",	 0, 0, 'v' }
		};

		c = getopt_long(argc, argv, "",
				long_options, &option_index);

		if(c == -1) {
			break;
		}

		switch(c) {
                case 'a':
                        options.alpha = atof(optarg);
                        break;
                case 'l':
                        options.lambda = atof(optarg);
                        break;
                case 'e':
                        options.save = string(optarg);
                        break;
                case 's':
                        options.samples = atoi(token(string(optarg), 0).c_str());
                        options.burnin  = atoi(token(string(optarg), 1).c_str());
                        break;
                case 't':
                        options.tfbs_length = atoi(optarg);
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
