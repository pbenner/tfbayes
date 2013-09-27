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
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <assert.h>

#include <fstream>

#include <tfbayes/alignment/alignment.hh>
#include <tfbayes/interface/datatypes.hh>
#include <tfbayes/dpm/init.hh>
#include <tfbayes/dpm/dpm-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-io.hh>
#include <tfbayes/dpm/dpm-tfbs-mean.hh>
#include <tfbayes/dpm/dpm-tfbs-sampler.hh>
#include <tfbayes/dpm/dpm-partition.hh>
#include <tfbayes/dpm/utility.hh>
#include <tfbayes/exception/exception.h>
#include <tfbayes/utility/linalg.hh>

using namespace std;

// options and global variables
// -----------------------------------------------------------------------------

typedef struct {
        size_t tfbs_length;
        double alpha;
        double discount;
        double lambda;
        double initial_temperature;
        const char* process_prior;
        const char* background_model;
        matrix_t* background_alpha;
        size_t background_context;
        const char* background_weights;
        vector_t*  baseline_weights;
        matrix_t** baseline_priors;
        const char** baseline_tags;
        size_t baseline_n;
        dpm_partition_t* partition;
        size_t population_size;
        const char* socket_file;
} options_t;

static options_t _options;
static dpm_tfbs_pmcmc_t* _sampler;
static sequence_data_t<data_tfbs_t::code_t> _phylogenetic_data;
static alignment_set_t<short> _alignment_set;

ostream&
operator<<(std::ostream& o, const options_t& options) {
        o << "Options:"              << endl
          << "-> alpha               = " << options.alpha               << endl
          << "-> discount            = " << options.discount            << endl
          << "-> lambda              = " << options.lambda              << endl
          << "-> initial temperature = " << options.initial_temperature << endl
          << "-> tfbs_length         = " << options.tfbs_length         << endl
          << "-> process prior       = " << options.process_prior       << endl
          << "-> background model    = " << options.background_model    << endl
          << "-> background context  = " << options.background_context  << endl
          << "-> background weights  = " << options.background_weights  << endl
          << "-> population_size     = " << options.population_size     << endl
          << "-> socket_file         = " << options.socket_file         << endl;
        return o;
}

__BEGIN_DECLS

// python interface
// -----------------------------------------------------------------------------

options_t* _dpm_tfbs_options()
{
        return &_options;
}

void _dpm_tfbs_init(const char* phylogenetic_input, const char* alignment_input)
{
        tfbs_options_t tfbs_options;

        __dpm_init__();

        _phylogenetic_data = data_tfbs_t::read_fasta(phylogenetic_input);
        _alignment_set = alignment_set_t<short>(alignment_input);

        // background alpha
        for (size_t i = 0; i < _options.background_alpha->rows; i++) {
                tfbs_options.background_alpha.push_back(
                        vector<double>(_options.background_alpha->columns, 0));
                for (size_t j = 0; j < _options.background_alpha->columns; j++) {
                        tfbs_options.background_alpha[i][j] = 
                                _options.background_alpha->mat[i][j];
                        cout << tfbs_options.background_alpha[i][j]
                             << endl;
                }
        }

        // baseline priors, names, and weights
        for (size_t k = 0; k < _options.baseline_n; k++) {
                tfbs_options.baseline_weights.push_back(_options.baseline_weights->vec[k]);
                tfbs_options.baseline_priors.push_back(matrix<double>());
                for (size_t i = 0; i < _options.baseline_priors[k]->rows; i++) {
                        tfbs_options.baseline_priors[k].push_back(
                                vector<double>(_options.baseline_priors[k]->columns, 0));
                        for (size_t j = 0; j < _options.baseline_priors[k]->columns; j++) {
                                tfbs_options.baseline_priors[k][i][j] = 
                                        _options.baseline_priors[k]->mat[i][j];
                        }
                }
                tfbs_options.baseline_tags.push_back(_options.baseline_tags[k]);
        }

        // tfbs options
        tfbs_options.alpha               = _options.alpha;
        tfbs_options.discount            = _options.discount;
        tfbs_options.lambda              = _options.lambda;
        tfbs_options.initial_temperature = _options.initial_temperature;
        tfbs_options.tfbs_length         = _options.tfbs_length;
        tfbs_options.process_prior       = _options.process_prior;
        tfbs_options.background_model    = _options.background_model;
        tfbs_options.background_context  = _options.background_context;
        tfbs_options.background_weights  = _options.background_weights;
        tfbs_options.socket_file         = _options.socket_file;

        if (_options.partition) {
                tfbs_options.partition   = *_options.partition;
        }

        _sampler = new dpm_tfbs_pmcmc_t(tfbs_options, _phylogenetic_data, _alignment_set, _options.population_size);

        cout << _options << endl;
}

void _dpm_tfbs_save(const char* filename)
{
        if (filename == NULL) {
                dpm_tfbs_save_result(cout, *_sampler);
        }
        else {
                ofstream file;
                file.open(filename);
                dpm_tfbs_save_result(file, *_sampler);
                file.close();
        }
}

unsigned int _dpm_tfbs_num_clusters() {
        return _sampler->gdpm()[0]->state().size();
}

matrix_t* _dpm_tfbs_get_posterior() {
        matrix_t* result;
        const vector<vector<double> >& probabilities = _sampler->samples().probabilities;
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
        size_t n = _sampler->data().size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < _sampler->data().size(i)) {
                        m = _sampler->data().size(i);
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
        for (data_tfbs_t::const_iterator it = _sampler->data().begin();
             it != _sampler->data().end(); it++) {
                const index_i& index = **it;
                result->mat[index[0]][index[1]] = _sampler->gdpm()[0]->state()[index];
        }
        return result;
}

vector_t* _dpm_tfbs_hist_likelihood() {
        const vector<double>& likelihood = _sampler->sampling_history().likelihood[0];
        size_t size = likelihood.size();
        vector_t* result = alloc_vector(size);

        for (size_t i = 0; i < size; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

vector_t* _dpm_tfbs_hist_switches() {
        const vector<double>& switches = _sampler->sampling_history().switches[0];
        size_t size = switches.size();
        vector_t* result = alloc_vector(size);

        for (size_t i = 0; i < size; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

void _dpm_tfbs_print() {
        cout << _sampler->gdpm()[0] << endl;
}

void _dpm_tfbs_sample(unsigned int n, unsigned int burnin) {
        _sampler->sample((size_t)n, (size_t)burnin);
}

void _dpm_tfbs_optimize() {
        _sampler->optimize();
}

void _dpm_tfbs_free() {
        delete(_sampler);
}

// handling partitions to resume a previous sampling state
// -----------------------------------------------------------------------------

dpm_partition_t* _dpm_partition_new()
{
        return new dpm_partition_t();
}

void _dpm_partition_add_component(dpm_partition_t* partition, const char* subset_tag)
{
        partition->add_component(subset_tag);
}

void _dpm_partition_add_index(dpm_partition_t* partition, int sequence, int position)
{
        seq_index_t index(sequence, position);

        partition->back().insert(index);
}

void _dpm_partition_free(dpm_partition_t* partition)
{
        delete(partition);
}

// handling lists of partitions to compute means and medians
// -----------------------------------------------------------------------------

vector<dpm_partition_t>* _dpm_partition_list_new()
{
        return new vector<dpm_partition_t>();
}

void _dpm_partition_list_add_partition(vector<dpm_partition_t>* partition_list, const dpm_partition_t* partition)
{
        partition_list->push_back(*partition);
}

void _dpm_partition_list_free(vector<dpm_partition_t>* partition_list)
{
        delete(partition_list);
}


// compute means and medians
// -----------------------------------------------------------------------------

dpm_partition_t* _dpm_mean(vector<dpm_partition_t>* partition_list)
{
        const dpm_partition_t& result = dpm_tfbs_mean(*partition_list, _phylogenetic_data.sizes(),
                                                      _options.tfbs_length, true);

        return new dpm_partition_t(result);
}

dpm_partition_t* _dpm_median(vector<dpm_partition_t>* partition_list)
{
        const dpm_partition_t& result = dpm_tfbs_median(*partition_list, _phylogenetic_data.sizes(),
                                                        _options.tfbs_length, true);

        return new dpm_partition_t(result);
}

__END_DECLS
