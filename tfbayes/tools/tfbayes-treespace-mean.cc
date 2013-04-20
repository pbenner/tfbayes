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
        bool cyclic;
        size_t drop;
        size_t k;
        size_t iterations;
        _options_t()
                : cyclic(false),
                  drop(0),
                  k(1),
                  iterations(100)
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
                      "Commands: mean, median\n"
                      "\n"
                      "Options:\n"
                      "             -c              - use cyclic version\n"
                      "             -n INTEGER      - number of iterations\n"
                      "             -d INTEGER      - drop first n trees\n"
                      "             -k INTEGER      - compute mean from every kth tree\n"
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
        list<pt_root_t*>  tree_list = parse_tree_list();
        list<ntree_t   > ntree_list;
        // convert trees
        list<pt_root_t*>::const_iterator it;
        size_t i;
        for (it = tree_list.begin(), i = 0; it != tree_list.end(); it++, i++) {
                pt_root_t* pt_root = *it;
                if (i >= options.drop && i % options.k == 0) {
                        ntree_list.push_back(pt_root);
                }
                pt_root->destroy();
        }

        return ntree_list;
}

void mean(const string& command)
{
        /* init random number generator */
        init();

        ntree_t result;
        /* phylogenetic tree */
        list<ntree_t> ntree_list = parse_tree_file();
        /* return if there is no tree in the list */
        if (ntree_list.size() == 0) return;

        if (command == "mean") {
                result = options.cyclic ?
                        mean_tree_cyc (ntree_list, options.iterations) :
                        mean_tree_rand(ntree_list, options.iterations);
        }
        else if (command == "median") {
                result = options.cyclic ?
                        median_tree_cyc (ntree_list, options.iterations) :
                        median_tree_rand(ntree_list, options.iterations);
        }
        else {
                wrong_usage("Unknown command.");
        }
        /* print result */
        pt_root_t* tmp = result.export_tree();
        cout << newick_format(tmp) << endl;
        tmp->destroy();
}

int main(int argc, char *argv[])
{
        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'v' }
                };

                c = getopt_long(argc, argv, "cd:k:n:",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'c':
                        options.cyclic = true;
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
        string command(argv[optind]);

        mean(command);

        return 0;
}
