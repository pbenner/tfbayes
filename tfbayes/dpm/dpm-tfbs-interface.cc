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
#include <pmcmc.hh>

#include <tfbayes/exception.h>
#include <tfbayes/fasta.hh>
#include <tfbayes/linalg.h>

using namespace std;

// options and global variables
// -----------------------------------------------------------------------------

typedef struct : public tfbs_options_t {
        size_t population_size;
} options_t;

static options_t _options;
static DpmTfbs* _gdpm;
static data_tfbs_t* _data;
static data_tfbs_t* _data_comp;
static GibbsSampler* _sampler;
static PopulationMCMC* _pmcmc;
static vector<string> _sequences;
static vector<string> _sequences_comp;

ostream&
operator<<(std::ostream& o, const options_t& options) {
        o << "Options:"              << endl
          << "-> tfbs_length     = " << options.tfbs_length     << endl
          << "-> process prior   = " << options.process_prior   << endl
          << "-> alpha           = " << options.alpha           << endl
          << "-> discount        = " << options.discount        << endl
          << "-> lambda          = " << options.lambda          << endl
          << "-> population_size = " << options.population_size << endl;
        return o;
}

// file i/o
// -----------------------------------------------------------------------------

ostream& operator<< (ostream& o, const ProductDirichlet& pd) {
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
        const posterior_t& posterior      = _pmcmc->posterior();
        const sampling_history_t& history = _pmcmc->sampling_history();
        const ClusterManager& cm          = _gdpm->clustermanager();

        file.setf(ios::showpoint);

        file << "[Result]" << endl;
        file << "posterior =" << endl;
        for (size_t i = 0; i < posterior.probabilities.size(); i++) {
                file << "\t";
                for (size_t j = 0; j < posterior.probabilities[i].size(); j++) {
                        file << (float)posterior.probabilities[i][j] << " ";
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
        for (Graph::const_iterator it = posterior.graph.begin();
             it != posterior.graph.end(); it++) {
                file << (*it).first.index1 << "-"
                     << (*it).first.index2 << "="
                     << static_cast<double>((*it).second)/static_cast<double>(_pmcmc->sampling_steps()) << " ";
        }
        file << endl;
        file << "hypergraph =" << endl;
        for (vector<string>::const_iterator it = posterior.hypergraph.begin();
             it != posterior.hypergraph.end(); it++) {
                file << "\t" << *it << endl;
        }
        for (ClusterManager::const_iterator it = cm.begin();
             it != cm.end(); it++) {
                if ((*it)->cluster_tag() == 0) {
                        file << "cluster_bg" << " =" << endl;
                        file << static_cast<const ProductDirichlet&>((*it)->model());
                }
        }
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
        _sequences_comp = complement(_sequences);

        // baseline priors
        std::vector<double> baseline_weights(_options.baseline_n, 0);
        gsl_matrix* baseline_priors[_options.baseline_n+1];
        for (size_t i = 0; i < _options.baseline_n; i++) {
                baseline_weights[i] = _options.baseline_weights->vec[i];
                baseline_priors[i]  = Simple::toGslMatrix(_options.baseline_priors[i]);
        }
        baseline_priors[_options.baseline_n] = NULL;

        _data      = new data_tfbs_t(_sequences, _options.tfbs_length);
        _gdpm      = new DpmTfbs(_options, *_data);
        _sampler   = new GibbsSampler(*_gdpm, *_data);
        _pmcmc     = new PopulationMCMC(*_sampler, _options.population_size);

        for (size_t i = 0; i < _options.baseline_n; i++) {
                gsl_matrix_free(baseline_priors[i]);
        }

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
        return _gdpm->clustermanager().size();
}

Simple::Matrix* _dpm_tfbs_get_posterior() {
        Simple::Matrix* result;
        const vector<vector<double> >& probabilities = _gdpm->posterior().probabilities;
        size_t n = probabilities.size();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < probabilities[i].size()) {
                        m = probabilities[i].size();
                }
        }

        // allocate matrix
        result = Simple::allocMatrix(n, m);
        // copy posterior
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

Simple::Matrix* _dpm_tfbs_cluster_assignments() {
        Simple::Matrix* result;
        size_t n = _gdpm->data().length();
        size_t m = 0;

        // compute maximum length
        for (size_t i = 0; i < n; i++) {
                if (m < _gdpm->data().length(i)) {
                        m = _gdpm->data().length(i);
                }
        }

        // allocate matrix
        result = Simple::allocMatrix(n, m);
        // default initialization to -1
        for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < m; j++) {
                        result->mat[i][j] = -1;
                }
        }
        // copy posterior
        for (data_tfbs_t::const_iterator it = _gdpm->data().begin();
             it != _gdpm->data().end(); it++) {
                const index_t& index = **it;
                result->mat[index[0]][index[1]] = _gdpm->clustermanager()[index];
        }
        return result;
}

Simple::Vector* _dpm_tfbs_hist_likelihood() {
        const vector<double>& likelihood = _pmcmc->sampling_history().likelihood[0];
        size_t length = likelihood.size();
        Simple::Vector* result = Simple::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = likelihood[i];
        }

        return result;
}

Simple::Vector* _dpm_tfbs_hist_switches() {
        const vector<double>& switches = _pmcmc->sampling_history().switches[0];
        size_t length = switches.size();
        Simple::Vector* result = Simple::allocVector(length);

        for (size_t i = 0; i < length; i++) {
                result->vec[i] = switches[i];
        }

        return result;
}

void _dpm_tfbs_print() {
        cout << *_gdpm << endl;
}

void _dpm_tfbs_sample(unsigned int n, unsigned int burnin) {
        _pmcmc->sample((size_t)n, (size_t)burnin);
}

void _dpm_tfbs_free() {
        delete(_data);
        delete(_gdpm);
        delete(_pmcmc);
}

__END_DECLS
