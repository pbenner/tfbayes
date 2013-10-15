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

#include <getopt.h>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/phylotree/phylotree.hh>
#include <tfbayes/phylotree/phylotree-parser.hh>
#include <tfbayes/phylotree/phylotree-sampler.hh>
#include <tfbayes/phylotree/phylotree-gradient.hh>
#include <tfbayes/phylotree/phylotree-gradient-ascent.hh>
#include <tfbayes/exception/exception.h>

#define alphabet_size 5
typedef float code_t;

using namespace std;

class posterior_values {
public:
        posterior_values(const pt_pmcmc_hastings_t<alphabet_size, code_t>& mh)
                : mh(mh) { }

        std::ostream& operator()(std::ostream& o) const {
                if (mh.log_posterior_history.begin() ==
                    mh.log_posterior_history.end()) {
                        return o;
                }
                // print header
                o << "x";
                for (size_t i = 0; i < mh.log_posterior_history.size(); i++) {
                        o << " y" << (i+1);
                }
                o << endl;
                // print log posterior
                size_t n = mh.log_posterior_history.begin()->size();
                for (size_t i = 0; i < n; i++) {
                        o << (i+1);
                        for (std::list<vector<double> >::const_iterator it = mh.log_posterior_history.begin();
                             it != mh.log_posterior_history.end(); it++) {
                                o << " " << fixed << it->operator[](i);
                        }
                        o << endl;
                }
                return o;
        }
protected:
        const pt_pmcmc_hastings_t<alphabet_size, code_t>& mh;
};

ostream& operator<< (ostream& o, const posterior_values& pv)
{
        return pv(o);
}

ostream& operator<< (ostream& o, const pt_pmcmc_hastings_t<alphabet_size, code_t>& mh)
{
        // print trees
        for (std::list<pt_root_t>::const_iterator it = mh.samples.begin();
             it != mh.samples.end(); it++) {
                o << newick_format(*it) << endl;
        }

        return o;
}

size_t hash_value(const exponent_t<alphabet_size, code_t>& exponent)
{
        size_t seed = 0;
        seed += (size_t)exponent[0] << 0;
        seed += (size_t)exponent[1] << 2;
        seed += (size_t)exponent[2] << 4;
        seed += (size_t)exponent[3] << 6;
        seed += (size_t)exponent[4] << 8;

        return seed;
}

void init() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

// Options
////////////////////////////////////////////////////////////////////////////////

typedef struct _options_t {
        double alpha;
        size_t burnin;
        double lambda;
        double r;
        size_t max_steps;
        double min_change;
        double epsilon;
        double sigma;
        string posterior;
        size_t jobs;
        _options_t()
                : alpha(0.2),
                  burnin(1000),
                  lambda(0.1),
                  r(1.0),
                  max_steps(1000),
                  min_change(0.0005),
                  epsilon(0.001),
                  sigma(0.01),
                  posterior(""),
                  jobs(1)
                { }
} options_t;

static options_t options;

ostream&
operator<<(std::ostream& o, const options_t& options) {
        o << "Options:"              << endl
          << "-> alpha                 = " << options.alpha               << endl
          << "-> burnin                = " << options.burnin              << endl
          << "-> gamma scale           = " << options.lambda              << endl
          << "-> gamma shape           = " << options.r                   << endl
          << "-> max steps             = " << options.max_steps           << endl
          << "-> min change            = " << options.min_change          << endl
          << "-> gradient step size    = " << options.epsilon             << endl
          << "-> proposal variance     = " << options.sigma               << endl
          << "-> save posterior values = " << options.posterior           << endl;
        return o;
}

// Main
////////////////////////////////////////////////////////////////////////////////

static
void print_usage(char *pname, FILE *fp)
{
        (void)fprintf(fp, "\nUsage: %s [OPTION] METHOD TREE ALIGNMENT\n\n", pname);
        (void)fprintf(fp,
                      "Methods: gradient-ascent, metropolis-hastings\n"
                      "\n"
                      "Options:\n"
                      "             -a DOUBLE       - pseudo count\n"
                      "             -b INTEGER      - number of burn-in samples for metropolis-hastings\n"
                      "             -m INTEGER      - maximum number of steps\n"
                      "             -n DOUBLE       - stop gradient ascent if change is smaller than this value\n"
                      "             -e DOUBLE       - gradient ascent step size\n"
                      "             -r DOUBLE       - r parameter (shape) for the gamma prior\n"
                      "             -l DOUBLE       - lambda parameter (scale) for the gamma prior\n"
                      "             -p FILE         - save the value of the log posterior\n"
                      "                               for each sample to FILE\n"
                      "             -s DOUBLE       - sigma^2 parameter for proposal distribution\n"
                      "             -j INTEGER      - number of parallel jobs\n"
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

pt_root_t parse_tree_file(const char* file_tree)
{
        FILE* yyin = fopen(file_tree, "r");
        if (yyin == NULL) {
                std_err(PERR, "Could not open phylogenetic tree");
        }
        list<pt_root_t> tree_list = parse_tree_list(yyin);
        assert(tree_list.size() == 1);
        pt_root_t pt_root = tree_list.front();
        fclose(yyin);

        return pt_root;
}

void run_optimization(const string& method, const char* file_tree, const char* file_alignment)
{
        init();

        /* phylogenetic tree */
        pt_root_t pt_root = parse_tree_file(file_tree);

        /* pseudo counts */
        exponent_t<alphabet_size, code_t> alpha;
        for (size_t i = 0; i < alphabet_size; i++) {
                alpha[i] = options.alpha;
        }

        /* alignment */
        alignment_t<code_t> alignment(file_alignment, pt_root);
        assert(alignment.length() > 0);

        if (method == "gradient-ascent") {
                pt_gradient_ascent_t<alphabet_size, code_t> pt_gradient_ascent(pt_root, alignment, alpha, options.r, options.lambda, options.epsilon);
                pt_gradient_ascent.run(options.max_steps, options.min_change);
        }
        else if (method == "metropolis-hastings") {
                normal_jump_t jump(options.sigma);
//                gamma_jump_t jump(1.6, 0.4);
                pt_metropolis_hastings_t<alphabet_size, code_t> pt_metropolis_hastings(pt_root, alignment, alpha, options.r, options.lambda, jump, 0.5);
                pt_pmcmc_hastings_t<alphabet_size, code_t> pmcmc(options.jobs, pt_metropolis_hastings);
                pmcmc.sample(options.max_steps, options.burnin);
                /* print posterior values to separate file */
                if (options.posterior != "") {
                        ofstream csv(options.posterior.c_str());
                        if (!csv.is_open()) {
                                cerr << "Unable to open file: "
                                     << options.posterior
                                     << endl;
                                exit(EXIT_FAILURE);
                        }
                        csv << posterior_values(pmcmc) << endl;
                        csv.close();
                }
                /* print tree samples */
                cout << pmcmc;
        }
        else {
                cerr << "Unknown optimization method: " << method
                     << endl;
                exit(EXIT_FAILURE);
        }
}

int main(int argc, char *argv[])
{
        const char* file_tree;
        const char* file_alignment;

        for(;;) {
                int c, option_index = 0;
                static struct option long_options[] = {
                        { "help",            0, 0, 'h' },
                        { "version",         0, 0, 'v' }
                };

                c = getopt_long(argc, argv, "a:b:e:m:n:r:l:p:s:j:hv",
                                long_options, &option_index);

                if(c == -1) {
                        break;
                }

                switch(c) {
                case 'a':
                        options.alpha = atof(optarg);
                        break;
                case 'b':
                        options.burnin = atoi(optarg);
                        break;
                case 'e':
                        options.epsilon = atof(optarg);
                        break;
                case 'm':
                        options.max_steps = atoi(optarg);
                        break;
                case 'n':
                        options.min_change = atof(optarg);
                        break;
                case 'r':
                        options.r = atof(optarg);
                        break;
                case 'l':
                        options.lambda = atof(optarg);
                        break;
                case 's':
                        options.sigma = atof(optarg);
                        break;
                case 'p':
                        options.posterior = string(optarg);
                        break;
                case 'j':
                        options.jobs = atoi(optarg);
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
        if(optind+3 != argc) {
                wrong_usage("Wrong number of arguments.");
                exit(EXIT_FAILURE);
        }

        string method(argv[optind]);
        file_tree      = argv[optind+1];
        file_alignment = argv[optind+2];

        // print options
        cerr << options << endl;

        run_optimization(method, file_tree, file_alignment);

        return 0;
}
