/* Copyright (C) 2013, 2014 Philipp Benner
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
#include <tfbayes/phylotree/parser.hh>
#include <tfbayes/phylotree/treespace.hh>
#include <tfbayes/utility/progress.hh>
#include <tfbayes/utility/random.hh>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/unordered/unordered_map.hpp>

#define alphabet_size 5

typedef boost::unordered_map<topology_t, std::list<ntree_t> > topology_map_t;

using namespace std;

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        double cut;
        double credibility_level;
        bool random;
        bool variance;
        string mean_file;
        size_t drop;
        size_t k;
        size_t iterations;
        double step_size;
        size_t verbose;
        _options_t()
                : cut(1e-8),
                  credibility_level(0.95),
                  random(true),
                  variance(false),
                  mean_file(""),
                  drop(0),
                  k(1),
                  iterations(100),
                  step_size(1.0),
                  verbose(0)
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
                      "Commands:\n"
                      "\n"
                      "      credibility            - radius of credibility region\n"
                      "      mean                   - Frechet mean\n"
                      "      median                 - geometric median\n"
                      "      majority-consensus     - majority rule consensus tree with\n"
                      "                               average branch lengths\n"
                      "      simple-mean            - compute the mean for each topology\n"
                      "      variance               - Frechet variance\n"
                      "\n"
                      "Options:\n"
                      "             -c float        - remove edges from resulting tree\n"
                      "                               if the length is shorter than FLOAT\n"
                      "                               (default: %e)\n"
                      "             -r              - use random version of the algorithm [default]\n"
                      "             -y              - use cyclic version of the algorithm\n"
                      "             -m file         - provide the Frechet mean for computing\n"
                      "                               the Frechet variance or credibility region\n"
                      "             -n integer      - number of iterations\n"
                      "             -d integer      - drop first n trees\n"
                      "             -k integer      - compute mean from every kth tree\n"
                      "             -s float        - step size parameter\n"
                      "   --credibility-level float - level for computing the credibility region\n"
                      "                               (default: 0.95)\n"
                      "\n"
                      "             -v              - set verbose level to one\n"
                      "   --verbose integer         - set verbose level (from 0 to 4)\n"
                      "\n"
                      "   --help                    - print help and exit\n"
                      "   --version                 - print version information and exit\n\n",
                      options.cut);
}

void wrong_usage(const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `tfbayes-treespace-estimate --help' for more information.\n");

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

list<ntree_t>
parse_tree_file(
        const string& filename, size_t drop, size_t k,
        boost::optional<const pt_root_t&> ref_tree = boost::optional<const pt_root_t&>())
{
        list<pt_root_t>  tree_list = parse_tree_list(filename, drop, k, ref_tree);
        list<ntree_t  > ntree_list;
        // convert trees
        for (list<pt_root_t>::const_iterator it = tree_list.begin();
             it != tree_list.end(); it++) {
                ntree_list.push_back(*it);
        }

        return ntree_list;
}

list<ntree_t> randomize_ntree_list(const list<ntree_t>& ntree_list, boost::random::mt19937& gen)
{
        vector<ntree_t> tmp(ntree_list.begin(), ntree_list.end());
        // random generator for shuffling the tree list
        boost::uniform_int<> uni_dist;
        boost::variate_generator<boost::random::mt19937&, boost::uniform_int<> > rng(gen, uni_dist);

        random_shuffle(tmp.begin(), tmp.end(), rng);

        return list<ntree_t>(tmp.begin(), tmp.end());
}

struct topology_counts_t {
        topology_counts_t(const topology_t& topology, const size_t n)
                : topology(topology), n(n) { }
        topology_t topology;
        size_t n;
};
struct by_counts {
        bool operator()(const topology_counts_t& a, const topology_counts_t& b) { 
                return a.n > b.n;
        }
};

void simple_mean(list<ntree_t>& result, const list<ntree_t>& ntree_list)
{
        topology_map_t map;
        vector<topology_counts_t> counts;

        // construct a set of topologies
        for (list<ntree_t>::const_iterator it = ntree_list.begin();
             it != ntree_list.end(); it++) {
                map[it->topology()].push_back(*it);
        }
        for (topology_map_t::const_iterator it = map.begin();
             it != map.end(); it++) {
                counts.push_back(topology_counts_t(it->first, it->second.size()));
        }
        sort(counts.begin(), counts.end(), by_counts());
        // for each topology compute the mean
        for (size_t i = 0; i < counts.size(); i++) {
                if (options.verbose) {
                        cerr << progress_t(i/(double)counts.size());
                }
                result.push_back(mean_same_topology(map[counts[i].topology]));
        }
}

void estimate(const string& command)
{
        // random number generator
        ////////////////////////////////////////////////////////////////////////
        boost::random::mt19937 gen;
        /* init random number generator */
        seed_rng(gen);

        // if variance should be computed take the mean tree
        // as the reference tree
        ////////////////////////////////////////////////////////////////////////
        boost::optional<pt_root_t> ref_tree;
        if (command == "variance" || command == "credibility") {
                /* read Frechet mean from file */
                if (options.mean_file == "") {
                        wrong_usage("Please provide the Frechet mean.");
                }
                list<pt_root_t> tmp = parse_tree_list(options.mean_file, 0, 1);
                assert(tmp.size() == 1);
                ref_tree = tmp.front();
        }

        // parse tree list
        ////////////////////////////////////////////////////////////////////////
        list<ntree_t> result_list;
        list<ntree_t> ntree_list;
        /* phylogenetic tree */
        if (ref_tree) {
                ntree_list = randomize_ntree_list(
                        parse_tree_file("", options.drop, options.k, *ref_tree), gen);
        } else {
                ntree_list = randomize_ntree_list(
                        parse_tree_file("", options.drop, options.k), gen);
        }
        /* return if there is no tree in the list */
        if (ntree_list.size() == 0) return;

        // run command
        ////////////////////////////////////////////////////////////////////////
        if (command == "credibility") {
                cout << frechet_credibility(ntree_list, ntree_t(*ref_tree), options.credibility_level)
                     << endl;
        }
        else if (command == "mean") {
                result_list.push_back(
                        options.random ?
                        mean_tree_rand(ntree_list, options.iterations, gen, default_lambda_t(options.step_size), options.verbose) :
                        mean_tree_cyc (ntree_list, options.iterations, default_lambda_t(options.step_size), options.verbose));
        }
        else if (command == "median") {
                result_list.push_back(
                        options.random ?
                        median_tree_rand(ntree_list, options.iterations, gen, default_lambda_t(options.step_size), options.verbose) :
                        median_tree_cyc (ntree_list, options.iterations, default_lambda_t(options.step_size), options.verbose));
        }
        else if (command == "majority-consensus") {
                result_list.push_back(
                        majority_consensus(ntree_list, options.verbose));
        }
        else if (command == "simple-mean") {
                simple_mean(result_list, ntree_list);
        }
        else if (command == "variance") {
                cout << frechet_variance(ntree_list, ntree_t(*ref_tree))
                     << endl;
        }
        else {
                wrong_usage("Unknown command.");
        }
        /* remove edges that are too short */
        for (list<ntree_t>::iterator it = result_list.begin(); it != result_list.end(); it++) {
                for (nedge_set_t::iterator is = it->nedge_set().begin();
                     is != it->nedge_set().end(); is++) {
                        if (is->d() < options.cut) {
                                is = it->nedge_set().erase(is);
                        }
                        // if the last element is removed from the
                        // edge set, then this prevents further
                        // incrementation of the iterator
                        if (is == it->nedge_set().end()) {
                                break;
                        }
                }
        }
        /* print resulting tree */
        for (list<ntree_t>::iterator it = result_list.begin(); it != result_list.end(); it++) {
                pt_root_t tmp = it->export_tree();
                if (options.verbose) {
                        cerr << *it << endl;
                }
                cout << newick_format(tmp) << endl;
        }

        /* free glp library space */
        glp_free_env();
}

int main(int argc, char *argv[])
{
        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "credibility-level", 1, 0, 'a' },
                        { "verbose",           1, 0, 'v' },
                        { "help",              0, 0, 'h' },
                        { "version",           0, 0, 'q' }
                };

                c = getopt_long(argc, argv, "c:rym:d:k:n:s:v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        if (atof(optarg) < 0 || 1 < atof(optarg)) {
                                print_usage(argv[0], stdout);
                                exit(EXIT_SUCCESS);
                        }
                        options.credibility_level = atof(optarg);
                        break;
                case 'c':
                        options.cut = atof(optarg);
                        break;
                case 'r':
                        options.random = true;
                        break;
                case 'y':
                        options.random = false;
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
                case 'm':
                        options.mean_file = string(optarg);
                        break;
                case 'n':
                        options.iterations = atoi(optarg);
                        break;
                case 's':
                        if (atoi(optarg) < 1.0) {
                                print_usage(argv[0], stdout);
                                exit(EXIT_SUCCESS);
                        }
                        options.step_size = atof(optarg);
                        break;
                case 'v':
                        if (optarg) {
                                options.verbose = atoi(optarg);
                        }
                        else {
                                options.verbose = 1;
                        }
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

        estimate(command);

        return 0;
}
