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
#include <map>
#include <algorithm>
#include <vector>
#include <sys/time.h>

#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/treespace.hh>
#include <tfbayes/exception/exception.h>

#include <boost/unordered/unordered_map.hpp>

using namespace std;

#define alphabet_size 5
typedef float code_t;

typedef std::map<named_nsplit_t, std::vector<double> > split_map_t;
typedef boost::unordered_map<topology_t, size_t> topology_map_t;

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        size_t drop;
        size_t k;
        bool verbose;
        _options_t()
                : drop(0),
                  k(1),
                  verbose(false)
                { }
} options_t;

static options_t options;

// Main
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION]... COMMAND < TREE_LIST\n\n", pname);
        (void)fprintf(fp,
                      "Format the list of posterior samples such that a histogram\n"
                      "can be easily computed.\n"
                      "\n"
                      "Commands:\n"
                      "             edges           - histogram of edge lengths\n"
                      "             topology        - histogram of topologies\n"
                      "\n"
                      "Options:\n"
                      "             -d INTEGER      - drop first n trees\n"
                      "             -k INTEGER      - compute histogram from every kth tree\n"
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
                      "Try `tfbayes-treespace-histogram --help' for more information.\n");

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

ostream& operator<<(ostream& o, const split_map_t& map)
{
        size_t size = 0;
        // determine size
        for (split_map_t::const_iterator it = map.begin(); it != map.end(); it++) {
                if (it->second.size() > size) {
                        size = it->second.size();
                }
        }
        // print splits
        size_t i = 0;
        for (split_map_t::const_iterator it = map.begin(); it != map.end(); it++, i++) {
                cout << "# s" << i << ": "
                     << it->first
                     << endl;
        }
        // print header
        for (size_t i = 0; i < map.size(); i++) {
                if (i != 0) o << " ";
                o << "s" << i;
        }
        o << endl;
        // print table
        for (size_t i = 0; i < size; i++) {
                for (split_map_t::const_iterator it = map.begin(); it != map.end(); it++) {
                        if (it != map.begin()) o << " ";
                        if (it->second.size() >= i) {
                                o << fixed << it->second[i];
                        }
                        else {
                                o << "      NA";
                        }
                }
                o << endl;
        }
        return o;
}

ostream& operator<<(ostream& o, const topology_map_t& map)
{
        vector<size_t> vec;

        for (topology_map_t::const_iterator it = map.begin();
             it != map.end(); it++) {
                vec.push_back(it->second);
        }
        sort(vec.begin(), vec.end(), greater<int>());

        o << "topology" << endl;
        for (size_t i = 0; i < vec.size(); i++) {
                for (size_t j = 0; j < vec[i]; j++) {
                        o << i+1 << endl;
                }
        }
        return o;
}

void histogram_edges()
{
        split_map_t map;
        ntree_t result;
        /* phylogenetic tree */
        list<ntree_t> ntree_list = parse_tree_file();
        /* return if there is no tree in the list */
        if (ntree_list.size() == 0) return;

        for (list<ntree_t>::const_iterator it = ntree_list.begin();
             it != ntree_list.end(); it++) {
                for (nedge_set_t::const_iterator is = it->nedge_set().begin();
                     is != it->nedge_set().end(); is++) {
                        const named_nsplit_t named_nsplit(**is, ntree_list.begin()->leaf_names());
                        map[named_nsplit].push_back(is->d());
                }
        }
        cout << map << endl;
}

void histogram_topology()
{
        topology_map_t map;

        ntree_t result;
        /* phylogenetic tree */
        list<ntree_t> ntree_list = parse_tree_file();
        /* return if there is no tree in the list */
        if (ntree_list.size() == 0) return;

        for (list<ntree_t>::const_iterator it = ntree_list.begin();
             it != ntree_list.end(); it++) {
                map[it->topology()]++;
        }
        cout << map << endl;
}

int main(int argc, char *argv[])
{
        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'q' }
                };

                c = getopt_long(argc, argv, "d:k:v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
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
        if (string(argv[optind]) == "edges") {
                histogram_edges();
        }
        else if (string(argv[optind]) == "topology") {
                histogram_topology();
        }
        else {
                wrong_usage("Unknown command.");
                exit(EXIT_FAILURE);
        }
        return 0;
}
