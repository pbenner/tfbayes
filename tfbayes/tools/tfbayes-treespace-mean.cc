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

#include <cstdio>
#include <getopt.h>
#include <sys/time.h>
#include <glpk.h>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/treespace.hh>
#include <tfbayes/exception/exception.h>

#define alphabet_size 5
typedef float code_t;

using namespace std;

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        bool random;
        bool variance;
        size_t drop;
        size_t k;
        size_t iterations;
        bool verbose;
        _options_t()
                : random(false),
                  variance(false),
                  drop(0),
                  k(1),
                  iterations(100),
                  verbose(false)
                { }
} options_t;

static options_t options;

// Main
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] COMMAND < TREE_LIST\n\n", pname);
        (void)fprintf(fp,
                      "Commands: mean, median, majority-consensus\n"
                      "\n"
                      "      mean                   - Frechet mean\n"
                      "      median                 - geometric median\n"
                      "      majority-consensus     - majority rule consensus tree with\n"
                      "                               average branch lengths\n"
                      "\n"
                      "Options:\n"
                      "             -r              - use random instead of cyclic version\n"
                      "             -f              - compute Frechet variance\n"
                      "             -n INTEGER      - number of iterations\n"
                      "             -d INTEGER      - drop first n trees\n"
                      "             -k INTEGER      - compute mean from every kth tree\n"
                      "             -v              - be verbose and print progress bar"
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
                      "Try `tfbayes-treespace-mean --help' for more information.\n");

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

list<ntree_t> parse_tree_file()
{
        list<pt_root_t*>  tree_list = parse_tree_list(NULL, options.drop, options.k);
        list<ntree_t   > ntree_list;
        // convert trees
        for (list<pt_root_t*>::const_iterator it = tree_list.begin();
             it != tree_list.end(); it++) {
                pt_root_t* pt_root = *it;
                ntree_list.push_back(pt_root);
                pt_root->destroy();
        }

        return ntree_list;
}

list<ntree_t> randomize_ntree_list(const list<ntree_t>& ntree_list)
{
        vector<ntree_t> tmp(ntree_list.begin(), ntree_list.end());

        random_shuffle(tmp.begin(), tmp.end());

        return list<ntree_t>(tmp.begin(), tmp.end());
}

void mean(const string& command)
{
        /* init random number generator */
        init();

        ntree_t result;
        /* phylogenetic tree */
        list<ntree_t> ntree_list = randomize_ntree_list(parse_tree_file());
        /* return if there is no tree in the list */
        if (ntree_list.size() == 0) return;

        if (command == "mean") {
                result = options.random ?
                        mean_tree_rand(ntree_list, options.iterations, default_lambda_t(), options.verbose) :
                        mean_tree_cyc (ntree_list, options.iterations, default_lambda_t(), options.verbose);
        }
        else if (command == "median") {
                result = options.random ?
                        median_tree_rand(ntree_list, options.iterations, default_lambda_t(), options.verbose) :
                        median_tree_cyc (ntree_list, options.iterations, default_lambda_t(), options.verbose);
        }
        else if (command == "majority-consensus") {
                result = majority_consensus(ntree_list, options.verbose);
        }
        else {
                wrong_usage("Unknown command.");
        }
        /* print resulting tree */
        pt_root_t* tmp = result.export_tree();
        cerr << result << endl;
        cout << newick_format(tmp) << endl;
        tmp->destroy();

        /* Frechet variance */
        if (command == "mean" && options.variance) {
                cerr << "Frechet variance: "
                     << frechet_variance(ntree_list, result)
                     << endl;
        }

        /* free glp library space */
        glp_free_env();
}

int main(int argc, char *argv[])
{
        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'q' }
                };

                c = getopt_long(argc, argv, "rfd:k:n:v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'r':
                        options.random = true;
                        break;
                case 'f':
                        options.variance = true;
                        break;
                case 'd':
                        if (atoi(optarg) < 0) {
                                print_usage(argv[0], stdout);
                                exit(EXIT_SUCCESS);
                        }
                        options.drop = atoi(optarg);
                        break;
                case 'k':
                        if (atoi(optarg) < 1) {
                                print_usage(argv[0], stdout);
                                exit(EXIT_SUCCESS);
                        }
                        options.k = atoi(optarg);
                        break;
                case 'n':
                        options.iterations = atoi(optarg);
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
        string command(argv[optind]);

        mean(command);

        return 0;
}
