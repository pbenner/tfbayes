/* Copyright (C) 2012, 2013 Philipp Benner
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
#include <fstream>
#include <cstdlib>
#include <ctime>

#include <getopt.h>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/fastarithmetics/fast-lnbeta.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-sampler.hh>
#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/phylotree/phylotree-gradient-ascent.hh>
#include <tfbayes/utility/strtools.hh>

#define alphabet_size 5
// type of the metropolis hastings algorithm
typedef pt_metropolis_hastings_t<alphabet_size> pt_mc_t;

using namespace std;

// format posterior values
////////////////////////////////////////////////////////////////////////////////

class posterior_values {
public:
        posterior_values(const pt_pmcmc_t<pt_mc_t>& mh)
                : mh(mh) { }

        std::ostream& operator()(std::ostream& o) const {
                if (mh.log_posterior_history().begin() ==
                    mh.log_posterior_history().end()) {
                        return o;
                }
                // print header
                o << "x";
                for (size_t i = 0; i < mh.log_posterior_history().size(); i++) {
                        o << " y" << (i+1);
                }
                o << endl;
                // print log posterior
                size_t n = mh.log_posterior_history().begin()->size();
                for (size_t i = 0; i < n; i++) {
                        o << (i+1);
                        for (std::list<vector<double> >::const_iterator it = mh.log_posterior_history().begin();
                             it != mh.log_posterior_history().end(); it++) {
                                o << " " << fixed << it->operator[](i);
                        }
                        o << endl;
                }
                return o;
        }
protected:
        const pt_pmcmc_t<pt_mc_t>& mh;
};

ostream& operator<< (ostream& o, const posterior_values& pv)
{
        return pv(o);
}

ostream& operator<< (ostream& o, const pt_pmcmc_t<pt_mc_t>& mh)
{
        // print trees
        for (std::list<pt_root_t>::const_iterator it = mh.samples().begin();
             it != mh.samples().end(); it++) {
                o << newick_format(*it) << endl;
        }

        return o;
}

// hashing functions
////////////////////////////////////////////////////////////////////////////////

size_t hash_value(const exponent_t<alphabet_size>& exponent)
{
        size_t seed = 0;
        seed += (size_t)exponent[0] << 0;
        seed += (size_t)exponent[1] << 2;
        seed += (size_t)exponent[2] << 4;
        seed += (size_t)exponent[3] << 6;
        seed += (size_t)exponent[4] << 8;

        return seed;
}

// options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        vector<double> alpha;
        double scale;
        double shape;
        size_t max_steps;
        double min_change;
        double step_size;
        double proposal_variance;
        string save_posterior;
        size_t chains;
        size_t threads;
        vector<double> temperatures;
        bool   verbose;
        _options_t()
                : alpha(alphabet_size, 0.2),
                  scale(0.1),
                  shape(1.0),
                  max_steps(1000),
                  min_change(0.0005),
                  step_size(0.001),
                  proposal_variance(0.01),
                  save_posterior(""),
                  chains(1),
                  threads(1),
                  temperatures(1,1),
                  verbose(false)
                { }
} options_t;

static options_t options;

ostream&
operator<<(std::ostream& o, const vector<double>& v)
{
        for (size_t i = 0; i < v.size(); i++) {
                if (i != 0) o << ":";
                o << v[i];
        }
        return o;
}

ostream&
operator<<(std::ostream& o, const options_t& options) {
        o << "Options:"                                                   << endl
          << "-> pseudocounts          = " << options.alpha               << endl
          << "-> gamma scale           = " << options.scale               << endl
          << "-> gamma shape           = " << options.shape               << endl
          << "-> max steps             = " << options.max_steps           << endl
          << "-> min change            = " << options.min_change          << endl
          << "-> gradient step size    = " << options.step_size           << endl
          << "-> proposal variance     = " << options.proposal_variance   << endl
          << "-> parallel chains       = " << options.chains              << endl
          << "-> number of threads     = " << options.threads             << endl
          << "-> temperatures          = " << options.temperatures        << endl
          << "-> save posterior values = " << options.save_posterior      << endl
          << "-> verbose               = " << options.verbose             << endl;
        return o;
}

// usage and version information
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] METHOD TREE ALIGNMENT\n\n", pname);
        (void)fprintf(fp,
                      "Methods: gradient-ascent, metropolis-hastings\n"
                      "\n"
                      "Options:\n"
                      "      --counts=f:f:f:f:f        - pseudo counts\n"
                      "      --steps=integer           - maximum number of steps\n"
                      "      --epsilon=float           - stop gradient ascent if change is smaller than this value\n"
                      "      --step-size=float         - gradient ascent step size\n"
                      "      --shape=float             - shape parameter for the gamma prior\n"
                      "      --scale=float             - scale parameter for the gamma prior\n"
                      "                                  for each sample to FILE\n"
                      "      --proposal-variance=float - variance of the proposal distribution\n"
                      "      --chains=integer          - number of parallel chains\n"
                      "      --threads=integer         - number of threads\n"
                      "      --save-posterior=file     - save the value of the log posterior\n"
                      "   -v                           - be verbose\n"
                      "\n"
                      "      --help                    - print help and exit\n"
                      "      --version                 - print version information and exit\n\n");
}

static
void wrong_usage(const char *msg)
{

        if(msg != NULL) {
                (void)fprintf(stderr, "%s\n", msg);
        }
        (void)fprintf(stderr,
                      "Try `tfbayes-treespace-optimize --help' for more information.\n");

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

static
pt_root_t parse_tree_file(const string& filename)
{
        list<pt_root_t> tree_list = parse_tree_list(filename);
        assert(tree_list.size() == 1);

        return tree_list.front();
}

// gradient ascent
////////////////////////////////////////////////////////////////////////////////

void run_gradient_ascent(
        const pt_root_t& pt_root,
        const alignment_t<>& alignment)
{
        pt_gradient_ascent_t<alphabet_size> pt_gradient_ascent(pt_root, alignment, options.alpha, options.shape, options.scale, options.step_size);
        pt_gradient_ascent.run(options.max_steps, options.min_change);
}

// sampler
////////////////////////////////////////////////////////////////////////////////

void save_posterior_values(const pt_pmcmc_t<pt_mc_t>& pmcmc)
{
        if (options.save_posterior != "") {
                ofstream csv(options.save_posterior.c_str());
                if (!csv.is_open()) {
                        cerr << "Unable to open file: "
                             << options.save_posterior
                             << endl;
                        exit(EXIT_FAILURE);
                }
                csv << posterior_values(pmcmc) << endl;
                csv.close();
        }
}

void run_mcmc(
        const pt_root_t& pt_root,
        const alignment_t<>& alignment)
{
        // a pool of threads for computing likelihoods
        thread_pool_t thread_pool(options.threads);
        // prior distribution on branch lengths
        boost::math::gamma_distribution<> gamma_distribution(options.shape, options.scale);
        // jumping distribution
        normal_jump_t jump(options.proposal_variance);
        // sampler
        pt_mc_t pt_mc(pt_root, alignment, options.alpha, gamma_distribution, jump, thread_pool);
        pt_pmcmc_t<pt_mc_t> pmcmc(options.chains, pt_mc);
        pmcmc(options.max_steps, options.verbose);
        // print posterior values to separate file
        save_posterior_values(pmcmc);
        // print tree samples
        cout << pmcmc;
}

void run_optimization(const string& method, const char* file_tree, const char* file_alignment)
{
        // phylogenetic tree
        pt_root_t pt_root = parse_tree_file(file_tree);

        // alignment
        alignment_t<> alignment(file_alignment, pt_root, nucleotide_alphabet_t(), options.verbose);
        assert(alignment.length() > 0);

        if (method == "gradient-ascent") {
                run_gradient_ascent(pt_root, alignment);
        }
        else if (method == "metropolis-hastings") {
                run_mcmc(pt_root, alignment);
        }
        else {
                cerr << "Unknown optimization method: " << method
                     << endl;
                exit(EXIT_FAILURE);
        }
}

// main function
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
        vector<string> tokens;
        const char* file_tree;
        const char* file_alignment;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "counts",            1, 0, 'a' },
                        { "chains",            1, 0, 'j' },
                        { "epsilon",           1, 0, 'n' },
                        { "proposal-variance", 1, 0, 's' },
                        { "shape",             1, 0, 'r' },
                        { "scale",             1, 0, 'l' },
                        { "steps",             1, 0, 'm' },
                        { "step-size",         1, 0, 'e' },
                        { "threads",           1, 0, 't' },
                        { "save-posterior",    1, 0, 'p' },
                        { "help",              0, 0, 'h' },
                        { "version",           0, 0, 'x' }
                };

                c = getopt_long(argc, argv, "v",
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
                case 'e':
                        options.step_size = atof(optarg);
                        break;
                case 'm':
                        options.max_steps = atoi(optarg);
                        break;
                case 'n':
                        options.min_change = atof(optarg);
                        break;
                case 'r':
                        options.shape = atof(optarg);
                        break;
                case 'l':
                        options.scale = atof(optarg);
                        break;
                case 's':
                        options.proposal_variance = atof(optarg);
                        break;
                case 'p':
                        options.save_posterior = string(optarg);
                        break;
                case 'j':
                        options.chains = atoi(optarg);
                        break;
                case 't':
                        options.threads = atoi(optarg);
                        break;
                case 'v':
                        options.verbose = true;
                        break;
                case 'h':
                        print_usage(argv[0], stdout);
                        exit(EXIT_SUCCESS);
                case 'x':
                        print_version(stdout);
                        exit(EXIT_SUCCESS);
                default:
                        wrong_usage(NULL);
                        exit(EXIT_FAILURE);
                }
        }
        if(optind+3 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }

        string method(argv[optind]);
        file_tree      = argv[optind+1];
        file_alignment = argv[optind+2];

        // print options
        if (options.verbose)
                cerr << options << endl;

        run_optimization(method, file_tree, file_alignment);

        return 0;
}
