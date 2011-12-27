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

#include <gsl/gsl_matrix.h>

#include <init.hh>
#include <dpm-tfbs-interface.hh>
#include <dpm-tfbs.hh>
#include <dpm-tfbs-sampler.hh>

#include <tfbayes/exception.h>
#include <tfbayes/fasta.hh>
#include <tfbayes/linalg.h>

using namespace std;

// options and global variables
// -----------------------------------------------------------------------------

typedef struct {
        size_t tfbs_length;
        double alpha;
        double discount;
        double lambda;
        size_t context;
        bool metropolis_optimize;
        const char* process_prior;
        const char* background_model;
        vector_t*  baseline_weights;
        matrix_t** baseline_priors;
        size_t baseline_n;
        size_t population_size;
} options_t;

static options_t _options;
static dpm_tfbs_sampler_t* _pmcmc;
static vector<string> _sequences;

ostream&
operator<<(std::ostream& o, const options_t& options) {
        o << "Options:"              << endl
          << "-> tfbs_length         = " << options.tfbs_length         << endl
          << "-> process prior       = " << options.process_prior       << endl
          << "-> background model    = " << options.background_model    << endl
          << "-> alpha               = " << options.alpha               << endl
          << "-> discount            = " << options.discount            << endl
          << "-> lambda              = " << options.lambda              << endl
          << "-> context             = " << options.context             << endl
          << "-> metropolis_optimize = " << options.metropolis_optimize << endl
          << "-> population_size     = " << options.population_size     << endl;
        return o;
}

// file i/o
// -----------------------------------------------------------------------------

ostream& operator<< (ostream& o, const product_dirichlet_t& pd) {
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
void read_file(const char* file_name, vector<string>& sequences)
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
}

static
void save_result(ostream& file)
{
        const samples_t& samples          = _pmcmc->samples();
        const sampling_history_t& history = _pmcmc->sampling_history();

        file.setf(ios::showpoint);

        file << "[Result]" << endl;
        file << "posterior =" << endl;
        for (size_t i = 0; i < samples.probabilities.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < samples.probabilities[i].size(); j++) {
                        file << (float)samples.probabilities[i][j] << " ";
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
        for (Graph::const_iterator it = samples.graph.begin();
             it != samples.graph.end(); it++) {
                file << *static_cast<const seq_index_t*>(&(*it).first.index1) << "-"
                     << *static_cast<const seq_index_t*>(&(*it).first.index2) << "="
                     << static_cast<double>((*it).second)/static_cast<double>(_pmcmc->sampling_steps()) << " ";
        }
        file << endl;
}

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

// python interface
// -----------------------------------------------------------------------------

options_t* _dpm_tfbs_options()
{
        return &_options;
}

void _dpm_tfbs_init(const char* filename)
{
        __dpm_init__();

        // read sequences
        read_file(filename, _sequences);

        // baseline priors
        vector<double> baseline_weights;
        vector<matrix<double> > baseline_priors;
        for (size_t k = 0; k < _options.baseline_n; k++) {
                baseline_weights.push_back(_options.baseline_weights->vec[k]);
                baseline_priors.push_back(matrix<double>());
                for (size_t i = 0; i < _options.baseline_priors[k]->rows; i++) {
                        baseline_priors[k].push_back(
                                vector<double>(_options.baseline_priors[k]->columns, 0));
                        for (size_t j = 0; j < _options.baseline_priors[k]->columns; j++) {
                                baseline_priors[k][i][j] = 
                                        _options.baseline_priors[k]->mat[i][j];
                        }
                }
        }

        // tfbs options
        tfbs_options_t tfbs_options;
        tfbs_options.alpha               = _options.alpha;
        tfbs_options.lambda              = _options.lambda;
        tfbs_options.discount            = _options.discount;
        tfbs_options.context             = _options.context;
        tfbs_options.metropolis_optimize = _options.metropolis_optimize;
        tfbs_options.tfbs_length         = _options.tfbs_length;
        tfbs_options.process_prior       = _options.process_prior;
        tfbs_options.background_model    = _options.background_model;
        tfbs_options.baseline_weights    = baseline_weights;
        tfbs_options.baseline_priors     = baseline_priors;

        _pmcmc = new dpm_tfbs_sampler_t(tfbs_options, _sequences, _options.population_size);

        cout << _options << endl;
}

void _dpm_tfbs_save(const char* filename)
{
        if (filename == NULL) {
                save_result(cout);
        }
        else {
                ofstream file;
                file.open(filename);
                save_result(file);
                file.close();
        }
}

unsigned int _dpm_tfbs_num_clusters() {
        return _pmcmc->_gdpm[0]->state().size();
}

matrix_t* _dpm_tfbs_get_posterior() {
        matrix_t* result;
        const vector<vector<double> >& probabilities = _pmcmc->samples().probabilities;
        size_t n = probabilities.size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < probabilities[i].size()) {
                        m = probabilities[i].size();
                }
        }

        // allocate matrix
        result = alloc_matrix(n, m);
        // copy samples
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < m; j++) {
                        if (j < probabilities[i].size()) {
                                result->mat[i][j] = probabilities[i][j];
                        }
                        else {
                                result->mat[i][j] = 0;
                        }
                }
        }

        return result;
}

matrix_t* _dpm_tfbs_cluster_assignments() {
        matrix_t* result;
        size_t n = _pmcmc->_data.size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < _pmcmc->_data.size(i)) {
                        m = _pmcmc->_data.size(i);
                }
        }

        // allocate matrix
        result = alloc_matrix(n, m);
        // default initialization to -1
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < m; j++) {
                        result->mat[i][j] = -1;
                }
        }
        // copy samples
        for (data_tfbs_t::const_iterator it = _pmcmc->_data.begin();
             it != _pmcmc->_data.end(); it++) {
                const index_i& index = **it;
                result->mat[index[0]][index[1]] = _pmcmc->_gdpm[0]->state()[index];
        }
        return result;
}

vector_t* _dpm_tfbs_hist_likelihood() {
        const vector<double>& likelihood = _pmcmc->sampling_history().likelihood[0];
        size_t size = likelihood.size();
        vector_t* result = alloc_vector(size);

        for (size_t i = 0; i < size; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

vector_t* _dpm_tfbs_hist_switches() {
        const vector<double>& switches = _pmcmc->sampling_history().switches[0];
        size_t size = switches.size();
        vector_t* result = alloc_vector(size);

        for (size_t i = 0; i < size; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

void _dpm_tfbs_print() {
        cout << _pmcmc->_gdpm[0] << endl;
}

void _dpm_tfbs_sample(unsigned int n, unsigned int burnin) {
        _pmcmc->sample((size_t)n, (size_t)burnin);
}

void _dpm_tfbs_free() {
        delete(_pmcmc);
}

__END_DECLS
