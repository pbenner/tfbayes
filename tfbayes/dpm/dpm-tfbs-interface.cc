/* Copyright (C) 2011-2013 Philipp Benner
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

ostream&
operator<<(std::ostream& o, const tfbs_options_t& options) {
        o << "Options:"              << endl
          << "-> alpha               = " << options.alpha               << endl
          << "-> discount            = " << options.discount            << endl
          << "-> lambda              = " << options.lambda              << endl
          << "-> initial temperature = " << options.initial_temperature << endl
          << "-> tfbs_length         = " << options.tfbs_length         << endl
          << "-> process prior       = " << *options.process_prior      << endl
          << "-> background model    = " << *options.background_model   << endl
          << "-> background context  = " << options.background_context  << endl
          << "-> background weights  = " << *options.background_weights << endl
          << "-> population_size     = " << options.population_size     << endl
          << "-> socket_file         = " << *options.socket_file        << endl;
        return o;
}

__BEGIN_DECLS

// python interface
// -----------------------------------------------------------------------------

void _dpm_init()
{
        __dpm_init__();
}

dpm_tfbs_pmcmc_t* _dpm_tfbs_init(const tfbs_options_t* options, const sampling_history_t* history)
{
        // print options
        cout << options << endl;
        // return new sampler
        return new dpm_tfbs_pmcmc_t(*options, *history);
}

void _dpm_tfbs_save(dpm_tfbs_pmcmc_t* sampler, const char* filename)
{
        if (filename == NULL) {
                dpm_tfbs_save_result(cout, *sampler);
        }
        else {
                ofstream file;
                file.open(filename);
                dpm_tfbs_save_result(file, *sampler);
                file.close();
        }
}

const sampling_history_t* _dpm_tfbs_results(dpm_tfbs_pmcmc_t* sampler)
{
        return &sampler->sampling_history();
}

void _dpm_tfbs_sample(dpm_tfbs_pmcmc_t* sampler, size_t n, size_t burnin) {
        sampler->sample(n, burnin);
}

void _dpm_tfbs_free(dpm_tfbs_pmcmc_t* sampler) {
        delete(sampler);
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

// compute point estimates
// -----------------------------------------------------------------------------

dpm_partition_t* _dpm_map(dpm_tfbs_pmcmc_t* sampler)
{
        return new dpm_partition_t(sampler->map());
}

dpm_partition_t* _dpm_mean(dpm_tfbs_pmcmc_t* sampler)
{
        return new dpm_partition_t(sampler->mean());
}

dpm_partition_t* _dpm_median(dpm_tfbs_pmcmc_t* sampler)
{
        return new dpm_partition_t(sampler->median());
}

__END_DECLS

// python interface
// -----------------------------------------------------------------------------

#include <locale>
#include <cctype>

#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(dpm_tfbs)
{
        class_<sampling_history_t>("sampling_history_t")
                ;
}
