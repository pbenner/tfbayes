/* Copyright (C) 2012-2014 Philipp Benner
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
#include <tfbayes/phylotree/parser.hh>
#include <tfbayes/phylotree/sampler.hh>
#include <tfbayes/phylotree/gradient.hh>
#include <tfbayes/phylotree/gradient-ascent.hh>
#include <tfbayes/utility/linalg.hh>
#include <tfbayes/utility/random.hh>
#include <tfbayes/utility/strtools.hh>

using namespace std;

// format posterior values
////////////////////////////////////////////////////////////////////////////////

class print_posterior_values {
public:
        explicit print_posterior_values(const pt_pmcmc_t& mh)
                : f(boost::bind(&print_posterior_values::print_mh, this, boost::cref(mh), _1))
                { }
        template <size_t AS, typename AC, typename PC>
        explicit print_posterior_values(const pt_gradient_ascent_t<AS, AC, PC>& ga)
                : f(boost::bind(&print_posterior_values::print_ga<AS, AC, PC>, this, boost::cref(ga), _1))
                { }

        ostream& operator()(ostream& o) const {
                return f(o);
        }
protected:
        template <size_t AS, typename AC, typename PC>
        ostream& print_ga(const pt_gradient_ascent_t<AS, AC, PC>& ga, ostream& o) const {
                size_t line = 0;
                // print header
                o << "x y" << endl;
                for (typename pt_gradient_ascent_t<AS, AC, PC>::history_t::values_t::const_iterator it = ga.history().values.begin();
                     it != ga.history().values.end(); it++) {
                        o << ++line << " " << *it << endl;
                }
                return o;
        }
        ostream& print_mh(const pt_pmcmc_t& mh, ostream& o) const {
                if (mh.history(0).values.begin() ==
                    mh.history(0).values.end()) {
                        return o;
                }
                size_t line = 0;
                // print header
                o << "x";
                for (size_t i = 0; i < mh.size(); i++) {
                        o << " y" << (i+1);
                }
                o << endl;
                // print the list such that the order of samples is preserved
                vector<pt_sampler_t::history_t::values_t::const_iterator> it_vec(mh.size());
                for (size_t i = 0; i < mh.size(); i++) {
                        it_vec[i] = mh.history(i).values.begin();
                }
                while (it_vec[0] != mh.history(0).values.end())
                {
                        o << (++line);
                        for (size_t i = 0; i < mh.size(); i++) {
                                o << " " << fixed << *it_vec[i];
                                // advance iteration for sampler i
                                it_vec[i]++;
                        }
                        o << endl;
                }
                return o;
        }
        boost::function<ostream& (ostream& o)> f;
};

class print_posterior_samples {
public:
        explicit print_posterior_samples(const pt_pmcmc_t& mh)
                : f(boost::bind(&print_posterior_samples::print_mh, this, boost::cref(mh), _1))
                { }
        template <size_t AS, typename AC, typename PC>
        explicit print_posterior_samples(const pt_gradient_ascent_t<AS, AC, PC>& ga)
                : f(boost::bind(&print_posterior_samples::print_ga<AS, AC, PC>, this, boost::cref(ga), _1))
                { }

        ostream& operator()(ostream& o) const {
                return f(o);
        }

protected:
        template <size_t AS, typename AC, typename PC>
        ostream& print_ga(const pt_gradient_ascent_t<AS, AC, PC>& ga, ostream& o) const {
                for (typename pt_gradient_ascent_t<AS, AC, PC>::history_t::samples_t::const_iterator it = ga.history().samples.begin();
                     it != ga.history().samples.end(); it++) {
                        o << newick_format(*it) << endl;
                }
                return o;
        }
        ostream& print_mh(const pt_pmcmc_t& mh, ostream& o) const {
                // print the list such that the order of samples is preserved
                vector<pt_sampler_t::history_t::samples_t::const_iterator> it_vec(mh.size());
                for (size_t i = 0; i < mh.size(); i++) {
                        it_vec[i] = mh.history(i).samples.begin();
                }
                while (it_vec[0] != mh.history(0).samples.end())
                {
                        for (size_t i = 0; i < mh.size(); i++) {
                                o << newick_format(*it_vec[i]) << endl;
                                // advance iteration for sampler i
                                it_vec[i]++;
                        }
                }
                return o;
        }
        boost::function<ostream& (ostream& o)> f;
};

ostream& operator<< (ostream& o, const print_posterior_values& pv)
{
        return pv(o);
}

ostream& operator<< (ostream& o, const print_posterior_samples& ps)
{
        return ps(o);
}

// options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        matrix<double> alpha;
        alphabet_t alphabet;
        double scale;
        double shape;
        size_t max_steps;
        double min_change;
        double leapfrog_step_size;
        size_t leapfrog_steps;
        double momentum_refreshment;
        double step_size;
        double proposal_variance;
        string save_posterior;
        size_t chains;
        size_t threads;
        vector<double> temperatures;
        bool   verbose;
        _options_t()
                : alphabet(nucleotide_alphabet_t()),
                  scale(0.1),
                  shape(1.0),
                  max_steps(1000),
                  min_change(0.0005),
                  leapfrog_step_size(0.01),
                  leapfrog_steps(1),
                  momentum_refreshment(0.0),
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
operator<<(ostream& o, const vector<double>& v)
{
        for (size_t i = 0; i < v.size(); i++) {
                if (i != 0) o << ":";
                o << v[i];
        }
        return o;
}

ostream&
operator<<(ostream& o, const matrix<double>& m)
{
        for (size_t i = 0; i < m.size(); i++) {
                if (i != 0) o << ";";
                o << m[i];
        }
        return o;
}

ostream&
operator<<(ostream& o, const options_t& options) {
        o << "Options:"                                                    << endl
          << "-> pseudocounts          = " << options.alpha                << endl
          << "-> gamma scale           = " << options.scale                << endl
          << "-> gamma shape           = " << options.shape                << endl
          << "-> leapfrog step_size    = " << options.leapfrog_step_size   << endl
          << "-> leapfrog steps        = " << options.leapfrog_steps       << endl
          << "-> momentum refreshment  = " << options.momentum_refreshment << endl
          << "-> max steps             = " << options.max_steps            << endl
          << "-> min change            = " << options.min_change           << endl
          << "-> gradient step size    = " << options.step_size            << endl
          << "-> proposal variance     = " << options.proposal_variance    << endl
          << "-> parallel chains       = " << options.chains               << endl
          << "-> number of threads     = " << options.threads              << endl
          << "-> temperatures          = " << options.temperatures         << endl
          << "-> save posterior values = " << options.save_posterior       << endl
          << "-> verbose               = " << options.verbose              << endl;
        return o;
}

// usage and version information
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] METHOD TREE ALIGNMENT\n\n", pname);
        (void)fprintf(fp,
                      "Methods: gradient-ascent          (only for testing purposes)\n"
                      "         metropolis-hastings\n"
                      "         hamiltonian-mcmc\n"
                      "\n"
                      "Options:\n"
                      "      --alphabet=string         - nucleotide (5) or protein (21) [default: nucleotide]\n"
                      "      --chains=integer          - number of parallel chains\n"
                      "      --counts=f:f:f:f:f[;...]  - matrix of pseudocounts\n"
                      "      --epsilon=float           - stop gradient ascent if change is smaller than this value\n"
                      "      --leapfrog-step-size      - integration step size for the hamiltonian mcmc\n"
                      "      --leapfrog-steps          - number of integration steps for the hamiltonian mcmc\n"
                      "      --proposal-variance=float - variance of the proposal distribution\n"
                      "      --hamiltonian-alpha=float - partial momentum refreshment parameter\n"
                      "                                  [-1,+1], where values closer to one cause a new sample\n"
                      "                                  to be similar to the old one\n"
                      "      --shape=float             - shape parameter for the gamma prior\n"
                      "      --scale=float             - scale parameter for the gamma prior\n"
                      "                                  for each sample to FILE\n"
                      "      --steps=integer           - maximum number of steps\n"
                      "      --step-size=float         - gradient ascent step size\n"
                      "      --threads=integer         - number of threads\n"
                      "      --temperatures=f:f:...    - a list of temperatures for the mc3\n"
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
                      "Try `tfbayes-treespace-sampler --help' for more information.\n");

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

template<typename T>
void save_posterior_values(const T& t)
{
        if (options.save_posterior != "") {
                ofstream csv(options.save_posterior.c_str());
                if (!csv.is_open()) {
                        cerr << "Unable to open file: "
                             << options.save_posterior
                             << endl;
                        exit(EXIT_FAILURE);
                }
                csv << print_posterior_values(t) << endl;
                csv.close();
        }
}

// gradient ascent
////////////////////////////////////////////////////////////////////////////////

// this gradient ascent algorithm is only ment for testing derivative
// computations, a more powerful algorithm for computing MAP solutions
// was developed by Felsenstein (i.e. a coordinate ascent algorithm)

template <size_t AS>
void run_gradient_ascent(
        const pt_root_t& pt_root,
        const alignment_map_t<>& alignment_map,
        const boost::math::gamma_distribution<>& gamma_distribution,
        thread_pool_t& thread_pool)
{
        pt_gradient_ascent_t<AS> pt_gradient_ascent(
                pt_root, alignment_map, options.alpha, gamma_distribution,
                thread_pool, options.step_size);
        pt_gradient_ascent(options.max_steps, options.min_change);
        // print posterior values to separate file
        save_posterior_values(pt_gradient_ascent);
        // print tree samples
        cout << print_posterior_samples(pt_gradient_ascent);
}

// sampler
////////////////////////////////////////////////////////////////////////////////

template <size_t AS>
void run_metropolis_hastings(
        const pt_root_t& pt_root,
        const alignment_map_t<>& alignment_map,
        const boost::math::gamma_distribution<>& gamma_distribution,
        thread_pool_t& thread_pool)
{
        typedef pt_metropolis_hastings_t<AS> pt_mc_t;

        threaded_rng_t rng;
        seed_rng(rng);
        // proposal distribution
        normal_proposal_t proposal(options.proposal_variance);
        // the metropolis sampler
        pt_mc_t pt_mc(pt_root, alignment_map, options.alpha, gamma_distribution, proposal, thread_pool);
        // parallel chains with different temperatures
        pt_mc3_t<pt_mc_t> pt_mc3(options.temperatures, pt_mc);
        // run several mc3 chains in parallel
        pt_pmcmc_t pmcmc(options.chains, pt_mc3);
        // execute the sampler
        pmcmc(options.max_steps, rng, options.verbose);
        // print posterior values to separate file
        save_posterior_values(pmcmc);
        // print tree samples
        cout << print_posterior_samples(pmcmc);
}

template <size_t AS>
void run_hamiltonian_mcmc(
        const pt_root_t& pt_root,
        const alignment_map_t<>& alignment_map,
        const boost::math::gamma_distribution<>& gamma_distribution,
        thread_pool_t& thread_pool)
{
        typedef pt_hamiltonian_t<AS> pt_ham_t;

        threaded_rng_t rng;
        seed_rng(rng);
        // the metropolis sampler
        pt_ham_t pt_ham(pt_root, alignment_map, options.alpha, gamma_distribution,
                        options.leapfrog_step_size, options.leapfrog_steps,
                        options.momentum_refreshment, thread_pool);
        // parallel chains with different temperatures
        pt_mc3_t<pt_ham_t> pt_mc3(options.temperatures, pt_ham);
        // run several mc3 chains in parallel
        pt_pmcmc_t pmcmc(options.chains, pt_mc3);
        // execute the sampler
        pmcmc(options.max_steps, rng, options.verbose);
        // print posterior values to separate file
        save_posterior_values(pmcmc);
        // print tree samples
        cout << print_posterior_samples(pmcmc);
}

void run_optimization(const string& method, const char* file_tree, const char* file_alignment)
{
        // a pool of threads for computing likelihoods
        thread_pool_t thread_pool(options.threads);
        // phylogenetic tree
        pt_root_t pt_root = parse_tree_file(file_tree);
        // prior distribution on branch lengths
        boost::math::gamma_distribution<> gamma_distribution(options.shape, options.scale);

        // alignment
        alignment_set_t<> alignment_set(file_alignment, pt_root, options.alphabet, options.verbose);
        assert(alignment_set.size() > 0);
        assert(alignment_set[0].length() > 0);
        // convert alignment
        alignment_map_t<> alignment_map(alignment_set);

        if (method == "gradient-ascent") {
                switch (options.alphabet.size()) {
                case  5: run_gradient_ascent< 5>(pt_root, alignment_map, gamma_distribution, thread_pool); break;
                case 21: run_gradient_ascent<21>(pt_root, alignment_map, gamma_distribution, thread_pool); break;
                default: assert(0);
                }
        }
        else if (method == "metropolis-hastings") {
                switch (options.alphabet.size()) {
                case  5: run_metropolis_hastings< 5>(pt_root, alignment_map, gamma_distribution, thread_pool); break;
                case 21: run_metropolis_hastings<21>(pt_root, alignment_map, gamma_distribution, thread_pool); break;
                default: assert(0);
                }
        }
        else if (method == "hamiltonian-mcmc") {
                switch (options.alphabet.size()) {
                case  5: run_hamiltonian_mcmc< 5>(pt_root, alignment_map, gamma_distribution, thread_pool); break;
                case 21: run_hamiltonian_mcmc<21>(pt_root, alignment_map, gamma_distribution, thread_pool); break;
                default: assert(0);
                }
        }
        else {
                cerr << "Unknown optimization method: " << method
                     << endl;
                exit(EXIT_FAILURE);
        }
}

// main function
////////////////////////////////////////////////////////////////////////////////

extern int optind, opterr, optopt;

int main(int argc, char *argv[])
{
        vector<string> tokens1;
        vector<string> tokens2;
        const char* file_tree;
        const char* file_alignment;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "alphabet",             1, 0, 'y' },
                        { "counts",               1, 0, 'a' },
                        { "chains",               1, 0, 'j' },
                        { "epsilon",              1, 0, 'n' },
                        { "proposal-variance",    1, 0, 's' },
                        { "leapfrog-step-size",   1, 0, 'b' },
                        { "leapfrog-steps",       1, 0, 'c' },
                        { "momentum-refreshment", 1, 0, 'd' },
                        { "shape",                1, 0, 'r' },
                        { "scale",                1, 0, 'l' },
                        { "steps",                1, 0, 'm' },
                        { "step-size",            1, 0, 'e' },
                        { "temperatures",         1, 0, 'u' },
                        { "threads",              1, 0, 't' },
                        { "save-posterior",       1, 0, 'p' },
                        { "help",                 0, 0, 'h' },
                        { "version",              0, 0, 'x' },
                        { 0,                      0, 0,  0  }
                };

                c = getopt_long(argc, argv, "v",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        tokens1 = token(string(optarg), ';');
                        for (size_t i = 0; i < tokens1.size(); i++) {
                                tokens2 = token(tokens1[i], ':');
                                options.alpha.push_back(vector<double>(tokens2.size(), 0.0));
                                for (size_t j = 0; j < tokens2.size(); j++) {
                                        options.alpha[i][j] = atof(tokens2[j].c_str());
                                }
                        }
                        break;
                case 'e':
                        options.step_size = atof(optarg);
                        break;
                case 'm':
                        options.max_steps = atoi(optarg);
                        break;
                case 'b':
                        options.leapfrog_step_size = atof(optarg);
                        break;
                case 'c':
                        options.leapfrog_steps = atoi(optarg);
                        break;
                case 'd':
                        options.momentum_refreshment = atof(optarg);
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
                case 'u':
                        tokens1 = token(string(optarg), ':');
                        options.temperatures = vector<double>();
                        for (size_t i = 0; i < tokens1.size(); i++) {
                                options.temperatures.push_back(atof(tokens1[i].c_str()));
                        }
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
                case 'y':
                        if (string("nucleotide") == string(optarg)) {
                                options.alphabet = nucleotide_alphabet_t();
                        }
                        else if (string("protein") == string(optarg)) {
                                options.alphabet = protein_alphabet_t();
                        }
                        else {
                                print_usage(argv[0], stdout);
                                exit(EXIT_FAILURE);
                        }
                        break;
                default:
                        wrong_usage(NULL);
                        exit(EXIT_FAILURE);
                }
        }
        if(optind+3 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }
        // generate prior if there was none specified
        if (options.alpha.size() == 0) {
                options.alpha = matrix<double>(1, options.alphabet.size(), 0.2);
        }
        // check dimensions of the prior
        for (size_t i = 0; i < options.alpha.size(); i++) {
                if (options.alpha[i].size() != options.alphabet.size()) {
                        wrong_usage("Prior does not match alphabet size.");
                        exit(EXIT_FAILURE);
                }
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
