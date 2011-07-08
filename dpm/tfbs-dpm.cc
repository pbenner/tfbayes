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

#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
//#include <cstdlib>
#include <ctime>

//#include <string.h>
#include <getopt.h>
#include <tfbayes/exception.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include <dpm.hh>
#include <init.hh>

using namespace std;

static
void print_usage(char *pname, FILE *fp)
{

	(void)fprintf(fp,
                      "\nUsage: %s [OPTION]... FILE\n\n", pname);
	(void)fprintf(fp,
                      "Options:\n"
                      "   --help	print help and exit\n"
                      "   --version	print version information and exit\n\n");

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
                        printf("%s\n", sequences[i]);
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

        fprintf(stderr, "max length: %lu, lines: %lu\n", (unsigned long)max_len, (unsigned long)lines);

        readfile(file_name, sequences);

        DPM*  gdpm = new DPM(lines, sequences);

        gdpm->gibbs_sample(10, 0);
        const vector<vector<double> >& posterior = gdpm->get_posterior();
        for (size_t i = 0; i < posterior.size(); i++) {
                for (size_t j = 0; j < posterior[i].size(); j++) {
                        printf("%f ", (float)posterior[i][j]);
                }
                cout << endl;
        }

        free(gdpm);
        free_sequences(sequences);
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
			{ "help",	0, 0, 'h' },
			{ "version",	0, 0, 'v' }
		};

		c = getopt_long(argc, argv, "",
				long_options, &option_index);

		if(c == -1) {
			break;
		}

		switch(c) {
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

	file_name = argv[optind];

        run_dpm(file_name);

        return 0;
}
