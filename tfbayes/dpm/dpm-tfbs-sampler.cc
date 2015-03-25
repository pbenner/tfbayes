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

#include <cmath> /* abs, ceil */

#include <boost/function.hpp>
#include <boost/format.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <tfbayes/dpm/dpm-tfbs-sampler.hh>
#include <tfbayes/utility/logarithmetic.hh>
#include <tfbayes/utility/statistics.hh>

using namespace std;
using namespace boost::asio;
using namespace boost::asio::local;

// dpm_tfbs_sampler_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(
        const tfbs_options_t& options,
        const dpm_tfbs_t& dpm_tfbs,
        const data_tfbs_t& data,
        save_queue_t<string>& output_queue)
        : gibbs_sampler_t        (dpm_tfbs, data, "", options.verbose)
        , phylogenetic_data      (&data)
        , m_output_queue         (&output_queue)
        , m_t0                   (options.initial_temperature)
        , m_block_samples        (options.block_samples)
        , m_block_samples_period (options.block_samples_period)
        , m_metropolis_proposals (options.metropolis_proposals)
        , m_optimize             (options.optimize)
        , m_optimize_period      (options.optimize_period)
        , m_verbose              (options.verbose)
{
        assert(options.initial_temperature  >= 1.0);
        assert(options.block_samples_period >= 1);
        assert(options.optimize_period >= 1);
}

dpm_tfbs_sampler_t::dpm_tfbs_sampler_t(const dpm_tfbs_sampler_t& sampler)
        : gibbs_sampler_t        (sampler)
        , phylogenetic_data      (sampler.phylogenetic_data)
        , m_output_queue         (sampler.m_output_queue)
        , m_t0                   (sampler.m_t0)
        , m_block_samples        (sampler.m_block_samples)
        , m_block_samples_period (sampler.m_block_samples_period)
        , m_metropolis_proposals (sampler.m_metropolis_proposals)
        , m_optimize             (sampler.m_optimize)
        , m_optimize_period      (sampler.m_optimize_period)
        , m_verbose              (sampler.m_verbose)
{ }

dpm_tfbs_sampler_t::~dpm_tfbs_sampler_t()
{ }

void
swap(dpm_tfbs_sampler_t& first, dpm_tfbs_sampler_t& second) {
        swap(static_cast<gibbs_sampler_t&>(first), static_cast<gibbs_sampler_t&>(second));
        swap(first.phylogenetic_data,      second.phylogenetic_data);
        swap(first.m_command_queue,        second.m_command_queue);
        swap(first.m_output_queue,         second.m_output_queue);
        swap(first.m_t0,                   second.m_t0);
        swap(first.m_block_samples,        second.m_block_samples);
        swap(first.m_block_samples_period, second.m_block_samples_period);
        swap(first.m_metropolis_proposals, second.m_metropolis_proposals);
        swap(first.m_optimize,             second.m_optimize);
        swap(first.m_optimize_period,      second.m_optimize_period);
        swap(first.m_verbose,              second.m_verbose);
}

dpm_tfbs_sampler_t*
dpm_tfbs_sampler_t::clone() const {
        return new dpm_tfbs_sampler_t(*this);
}

dpm_tfbs_sampler_t&
dpm_tfbs_sampler_t::operator=(const sampler_t& sampler)
{
        dpm_tfbs_sampler_t tmp(static_cast<const dpm_tfbs_sampler_t&>(sampler));
        swap(*this, tmp);
        return *this;
}

save_queue_t<command_t*>&
dpm_tfbs_sampler_t::command_queue()
{
        return m_command_queue;
}

const save_queue_t<command_t*>&
dpm_tfbs_sampler_t::command_queue() const
{
        return m_command_queue;
}

dpm_tfbs_t&
dpm_tfbs_sampler_t::dpm()
{
        return *static_cast<dpm_tfbs_t*>(m_dpm);
}

const dpm_tfbs_t&
dpm_tfbs_sampler_t::dpm() const
{
        return *static_cast<const dpm_tfbs_t*>(m_dpm);
}

// Gibbs samples
////////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
void
print_weights(
        size_t n,
        double log_weights1[],
        double log_weights2[],
        cluster_tag_t cluster_tags1[],
        cluster_tag_t cluster_tags2[])
{
        double ref = -std::numeric_limits<double>::infinity();

        cerr << "log weights: " << endl;
        for (size_t i = 0; i < n; i++) {
                cerr << boost::format("-> cluster %d has weight: %f (%f)") % cluster_tags1[i] % 
                        log(exp(log_weights1[i] - log_weights2[n-1]) - exp(ref - log_weights2[n-1])) % log_weights1[i]
                     << endl;
                ref = log_weights1[i];
        }
        for (size_t i = 0; i < n; i++) {
                cerr << boost::format("-> cluster %d has weight (reversed): %f (%f)") % cluster_tags2[i] %
                        log(exp(log_weights2[i] - log_weights2[n-1]) - exp(ref - log_weights2[n-1])) % log_weights2[i]
                     << endl;
                ref = log_weights2[i];
        }
}
#endif

bool
dpm_tfbs_sampler_t::m_gibbs_sample(const index_t& index, double temp, bool optimize) {
        size_t length;
        if (!dpm().state().get_free_range(index, length)) {
                return false;
        }
        range_t range1(index, length, false);
        range_t range2(index, length, true );
        ////////////////////////////////////////////////////////////////////////
        // first release the element from its cluster
        cluster_tag_t old_cluster_tag = dpm().state()[index];
        state().remove(range1);
        ////////////////////////////////////////////////////////////////////////
        size_t components = dpm().mixture_components() + dpm().baseline_components();
        double log_weights1[components];
        double log_weights2[components];
        cluster_tag_t cluster_tags1[components];
        cluster_tag_t cluster_tags2[components];
        ////////////////////////////////////////////////////////////////////////
        // compute weights
        dpm().mixture_weights(range1, log_weights1, cluster_tags1, temp, true);
        dpm().mixture_weights(range2, log_weights2, cluster_tags2, temp, false, log_weights1[components-1]);

        std::pair<size_t, size_t> result;
        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        if (optimize) {
                result = select_max_component2(components, log_weights1, log_weights2);
        }
        else {
                result = select_component2(components, log_weights1, log_weights2, gen());
        }
#ifdef DEBUG
        if (old_cluster_tag != cluster_tags1[result.second]) {
                cerr << "Moving " << static_cast<const index_t&>(index)
                     << " from "  << old_cluster_tag
                     << " to "    << cluster_tags1[result.second]
                     << endl;
                cerr << "Weights: "
                     << print_weights(components,
                                      log_weights1, log_weights2,
                                      cluster_tags1, cluster_tags2)
                     << endl;
                cerr << print_alignment_pretty(dpm().alignment_set()[range1])
                     << endl;
        }
#endif
        ////////////////////////////////////////////////////////////////////////
        if (result.first == 1) {
                dpm().state().add(range1, cluster_tags1[result.second]);
                return old_cluster_tag != cluster_tags1[result.second];
        }
        else {
                dpm().state().add(range2, cluster_tags2[result.second]);
                return old_cluster_tag != cluster_tags2[result.second];
        }
}

size_t
dpm_tfbs_sampler_t::m_gibbs_sample(double temp, bool optimize) {
        size_t sum = 0;
        // the indexer needs to be constant since it is shared between
        // processes, so to shuffle the indices we first need to
        // obtain a copy
        vector<index_t> indices(m_indexer->sampling_begin(), m_indexer->sampling_end());
        random_shuffle(indices.begin(), indices.end());
        // now sample
        for (vector<index_t>::const_iterator it = indices.begin();
             it != indices.end(); it++) {
                if(m_gibbs_sample(*it, temp, optimize)) sum+=1;
        }
        return sum;
}

// Gibbs block samples
////////////////////////////////////////////////////////////////////////////////

void
dpm_tfbs_sampler_t::m_block_sample(cluster_tag_t cluster_tag, double temp, bool optimize)
{
        cluster_t& cluster = state()[cluster_tag];
        vector<range_t> range_set;
        cluster_tag_t old_cluster_tag = cluster.cluster_tag();

        ////////////////////////////////////////////////////////////////////////
        if (cluster.size() == 0) {
                return;
        }
        ////////////////////////////////////////////////////////////////////////
        // fill range_set
        for (cluster_t::iterator it = cluster.begin(); it != cluster.end(); it++)
        {
                range_set.push_back(*it);
        }
        ////////////////////////////////////////////////////////////////////////
        // release all elemente from the cluster
        for (vector<range_t>::iterator it = range_set.begin(); it != range_set.end(); it++)
        {
                dpm().state().remove(*it);
        }
        ////////////////////////////////////////////////////////////////////////
        // obtian the mixture probabilities
        size_t components = dpm().mixture_components() + dpm().baseline_components();
        double log_weights[components];
        cluster_tag_t cluster_tags[components];
        dpm().mixture_weights(range_set, log_weights, cluster_tags, temp, true);

        cluster_tag_t new_cluster_tag;
        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        if (optimize) {
                new_cluster_tag = cluster_tags[select_max_component(components, log_weights)];
        }
        else {
                new_cluster_tag = cluster_tags[select_component(components, log_weights, gen())];
        }

        ////////////////////////////////////////////////////////////////////////
        // print some information to stderr
        flockfile(stderr);
        if (state()[new_cluster_tag].size() != 0 && m_verbose >= 2) {
                cerr << boost::format("%s: cluster %d merged with cluster %d (%d + %d)")
                        % m_name % old_cluster_tag % new_cluster_tag
                        % state()[new_cluster_tag].size() % range_set.size()
                     << endl;
        }
        fflush(stderr);
        funlockfile(stderr);

        ////////////////////////////////////////////////////////////////////////
        // move all elements to the new cluster
        for (size_t i = 0; i < range_set.size(); i++) {
                dpm().state().add(range_set[i], new_cluster_tag);
        }
}

void
dpm_tfbs_sampler_t::m_block_sample(double temp, bool optimize)
{
        // since clusters are modified it is not possible to simply
        // loop through the list of clusters, we need to be a bit more
        // careful here!
        vector<cluster_tag_t> used_clusters;
        for (cm_iterator it = dpm().state().begin(); it != dpm().state().end(); it++) {
                cluster_t& cluster = **it;
                if (!dpm().state().is_background(cluster)) {
                        used_clusters.push_back(cluster.cluster_tag());
                }
        }

        // go through the list of used clusters and if they are still
        // used then generate a block sample
        for (vector<cluster_tag_t>::const_iterator it = used_clusters.begin(); it != used_clusters.end(); it++) {
                m_block_sample(*it, temp, optimize);
        }
}

// Metropolis-Hastings samples
////////////////////////////////////////////////////////////////////////////////

bool
dpm_tfbs_sampler_t::m_metropolis_proposal_move(cluster_t& cluster, stringstream& ss)
{
        boost::random::uniform_int_distribution<> dist   (0, 1);
        boost::random::uniform_int_distribution<> dist_bg(0, dpm().state().bg_cluster_tags.size()-1);
        boost::random::uniform_int_distribution<> steps  (
                1, ceil(double(cluster.model().id().length)/2));
        size_t n = steps(gen());
        // select a background component at random
        cluster_tag_t bg_cluster_tag = dpm().state().bg_cluster_tags[dist_bg(gen())];
        // save old cluster size
        size_t size = cluster.size();

        dpm().state().save(cluster.cluster_tag(), bg_cluster_tag);

        if (dist(gen()) == 0) {
                ss << boost::format("%s: cluster %d: move %d steps to the right accepted (%d -> %d)")
                        % m_name % cluster.cluster_tag() % n % size % cluster.size();
                return dpm().state().move_right(cluster, bg_cluster_tag, n);
        }
        else {
                ss << boost::format("%s: cluster %d: move %d steps to the left accepted (%d -> %d)")
                        % m_name % cluster.cluster_tag() % n % size % cluster.size();
                return dpm().state().move_left(cluster, bg_cluster_tag, n);
        }
}

bool
dpm_tfbs_sampler_t::m_metropolis_proposal_size(cluster_t& cluster, stringstream& ss)
{
        boost::random::uniform_int_distribution<> dist   (0, 1);
        boost::random::uniform_int_distribution<> dist_bg(0, dpm().state().bg_cluster_tags.size()-1);

        // save the current state
        cluster_tag_t bg_cluster_tag = dpm().state().bg_cluster_tags[dist_bg(gen())];
        dpm().state().save(cluster.cluster_tag(), bg_cluster_tag);

        // get the foreground model
        const product_dirichlet_t& model = static_cast<const product_dirichlet_t&>(
                cluster.model());

        // find the position of the current length
        vector<size_t>::const_iterator it = find(
                model.lengths().begin(),
                model.lengths().end(),
                cluster.model().id().length);

        if (it == model.lengths().begin())
                goto increase;
        else if (it+1 == model.lengths().end())
                goto decrease;
        else if (dist(gen()) == 0)
                goto increase;
        else
                goto decrease;

increase:
        it++;
        ss << boost::format("%s: cluster %d: increasing length to %d accepted")
                % m_name % cluster.cluster_tag() % *it;
        return dpm().state().set_length(cluster, bg_cluster_tag, *it);
decrease:
        it--;
        ss << boost::format("%s: cluster %d: decreasing length to %d accepted")
                % m_name % cluster.cluster_tag() % *it;
        return dpm().state().set_length(cluster, bg_cluster_tag, *it);
}

bool
dpm_tfbs_sampler_t::m_metropolis_sample(cluster_tag_t cluster_tag, double temp, bool optimize,
                                        boost::function<bool (cluster_t& cluster, stringstream& ss)> f) {
        cluster_t& cluster = state()[cluster_tag];
        double posterior_ref = dpm().posterior();
        double posterior_tmp;
        stringstream ss;

        /* allocate a uniform distribution on the unit inverval */
        boost::random::uniform_01<> dist;

        if (f(cluster, ss)) {
                posterior_tmp = dpm().posterior();

                /* posterior value is on log scale! */
                if (optimize) {
                        if (posterior_tmp > posterior_ref) {
                                goto accepted;
                        }
                }
                else {
                        if (dist(gen()) <= min(exp((posterior_tmp - posterior_ref)/temp), 1.0)) {
                                goto accepted;
                        }
                }
                dpm().state().restore();
        }
        return false;

accepted:
        if (m_verbose >= 2) {
                flockfile(stderr);
                cerr << ss.str() << endl;
                fflush(stderr);
                funlockfile(stderr);
        }
        return true;
}

void
dpm_tfbs_sampler_t::m_metropolis_sample(
        double temp, bool optimize,
        boost::function<bool (cluster_t& cluster, stringstream& ss)> f)
{
        // the cluster list ist altered by metropolis samplers, hence
        // it is necessary to work with the list of cluster tags
        vector<cluster_tag_t> cluster_tags;
        for (cl_iterator it = state().begin(); it != state().end(); it++) {
                if (!dpm().state().is_background(**it)) {
                        cluster_tags.push_back((*it)->cluster_tag());
                }
        }
        for (vector<cluster_tag_t>::iterator it = cluster_tags.begin(); it != cluster_tags.end(); it++) {
                m_metropolis_sample(*it, temp, optimize, f);
        }
}

void
dpm_tfbs_sampler_t::m_metropolis_sample(double temp, bool optimize) {
        // sample multiple times
        for (size_t i = 0; i < m_metropolis_proposals; i++) {
                m_metropolis_sample(
                        temp, optimize, boost::bind(&dpm_tfbs_sampler_t::m_metropolis_proposal_size, this, _1, _2));
                m_metropolis_sample(
                        temp, optimize, boost::bind(&dpm_tfbs_sampler_t::m_metropolis_proposal_move, this, _1, _2));
        }
}

// Main
////////////////////////////////////////////////////////////////////////////////

size_t
dpm_tfbs_sampler_t::m_sample(size_t i, size_t n, double temp, bool optimize) {
        // call the standard hybrid sampler that first produces a
        // Gibbs sample and afterwards make a Metropolis-Hastings step
        size_t s = m_gibbs_sample(temp, optimize);
        m_metropolis_sample(temp, optimize);
        // do a Gibbs block sampling step, i.e. go through all
        // clusters and try to merge them
        if ((m_block_samples && i % m_block_samples_period == 0) || optimize) {
                m_block_sample(temp, optimize);
        }
        // update clusters
        for (cl_iterator it = dpm().state().begin(); it != dpm().state().end(); it++) {
                (**it).update();
        }
        if (m_verbose >= 3) {
                flockfile(stderr);
                if (m_verbose >= 2) {
                        BOOST_FOREACH(cluster_tag_t& tag, dpm().state().bg_cluster_tags) {
                                cerr << m_name << ": "
                                     << dpm().state()[tag].model().print_counts()
                                     << endl;
                        }
                }
                fflush(stderr);
                funlockfile(stderr);
        }
        // we are done with sampling here, now process commands
        while (!m_command_queue.empty()) {
                stringstream ss;
                command_t* command = m_command_queue.front();
                ss << command->operator()(dpm().state(), *this);
                m_command_queue.pop();
                delete(command);
                m_output_queue->push(ss.str());
        }
        return s;
}

size_t
dpm_tfbs_sampler_t::m_sample(size_t i, size_t n, bool is_burnin) {
        // temperature for simulated annealing
        size_t result = 0;
        double temp = 1.0;
        if (is_burnin) {
                // geometric decline of the temperature
                temp = m_t0*pow((1.0/m_t0), (double)i/n);
        }
        if (m_verbose >= 1) {
                flockfile(stderr);
                cerr << m_name << ": "
                     << "temperature is " << temp << endl;
                fflush(stderr);
                funlockfile(stderr);
        }
        // save temperature
        m_sampling_history.temperature[0].push_back(temp);
        if (!is_burnin && m_optimize && i % m_optimize_period == 0) {
                double old_posterior;
                double new_posterior = dpm().posterior();
                do {
                        result += m_sample(i, n, temp, true);
                        old_posterior = new_posterior;
                        new_posterior = dpm().posterior();
                        if (m_verbose >= 1) {
                                flockfile(stderr);
                                cerr << boost::format("%s: Optimizing... new posterior value: %f (increment: %f)")
                                        % m_name % new_posterior % (new_posterior - old_posterior)
                                     << endl;
                                fflush(stderr);
                                funlockfile(stderr);
                        }
                        // During optimization, local gibbs moves are
                        // used which in some cases might be
                        // incongruent with the global posterior
                        // probabilities, caused by the approximation
                        // of the gamma functions. Hence, we cannot
                        // assume that the posterior probability
                        // always increases!
                } while (abs(old_posterior - new_posterior) > 1e-4);
        }
        else {
                result = m_sample(i, n, temp, false);
        }
        return result;
}

void
dpm_tfbs_sampler_t::m_update_sampling_history(size_t switches)
{
        matrix<double> cluster_sizes = m_sampling_history.cluster_sizes;
        gibbs_sampler_t::m_update_sampling_history(switches);

        vector<double> tmp;

        BOOST_FOREACH(const cluster_t* cluster, dpm().state()) {
                if (!dpm().state().is_background(*cluster)) {
                        tmp.push_back(cluster->size());
                }
        }
        m_sampling_history.cluster_sizes = cluster_sizes;
        m_sampling_history.cluster_sizes.push_back(tmp);
}

// dpm_tfbs_pmcmc_t
////////////////////////////////////////////////////////////////////////////////

dpm_tfbs_pmcmc_t::dpm_tfbs_pmcmc_t(
        const tfbs_options_t& options,
        const sampling_history_t& history)
        : population_mcmc_t(options.population_size, history)
        , m_options        (options)
        , m_data           (options.phylogenetic_file)
        , m_alignment_set  (options.alignment_file, boost::optional<const pt_root_t&>(),
                            nucleotide_alphabet_t(), options.verbose)
        , m_socket_file    (options.socket_file)
        , m_server         (NULL)
        , m_bt             (NULL)
{
        // initialize dpm_tfbs
        dpm_tfbs_t dpm_tfbs(options, m_data, m_alignment_set);

        for (size_t i = 0; i < m_size; i++) {
                // check if we have to resume from an old state
                if (history.partitions.size() != 0) {
                        if (history.partitions.size() % m_size != 0) {
                                cerr << "Cannot resume from previous sampling run: Number of"
                                     << " partitions does not match."
                                     << endl;
                                exit(EXIT_FAILURE);
                        }
                        // resume partition
                        const dpm_partition_t& partition = history.partitions
                                [history.partitions.size()-m_size+i];
                        dpm_tfbs.state().set_partition(partition);
                }
                // initialize sampler
                m_population[i] = new dpm_tfbs_sampler_t(options, dpm_tfbs, m_data, m_output_queue);
                // make some noise
                std::stringstream ss;
                ss << "Sampler " << i+1;
                operator[](i).name() = ss.str();
        }
        m_start_server();
}

dpm_tfbs_pmcmc_t::dpm_tfbs_pmcmc_t(const dpm_tfbs_pmcmc_t& sampler)
        : population_mcmc_t (sampler)
        , m_data            (sampler.m_data)
{
        cerr << "Cannot copy dpm_tfbs_pmcmc_t!"
             << endl;
        exit(EXIT_FAILURE);
}

dpm_tfbs_pmcmc_t::~dpm_tfbs_pmcmc_t() {
        m_stop_server();
}

dpm_tfbs_pmcmc_t*
dpm_tfbs_pmcmc_t::clone() const {
        return new dpm_tfbs_pmcmc_t(*this);
}

void
swap(dpm_tfbs_pmcmc_t& first, dpm_tfbs_pmcmc_t& second) {
        swap(static_cast<population_mcmc_t&>(first),
             static_cast<population_mcmc_t&>(second));
        swap(first.m_options,       second.m_options);
        swap(first.m_data,          second.m_data);
        swap(first.m_alignment_set, second.m_alignment_set);
        swap(first.m_socket_file,   second.m_socket_file);
        swap(first.m_server,        second.m_server);
        swap(first.m_bt,            second.m_bt);
        swap(first.m_output_queue,  second.m_output_queue);
}

dpm_tfbs_pmcmc_t&
dpm_tfbs_pmcmc_t::operator=(const sampler_t& sampler)
{
        dpm_tfbs_pmcmc_t tmp(static_cast<const dpm_tfbs_pmcmc_t&>(sampler));
        swap(*this, tmp);
        return *this;
}

const tfbs_options_t&
dpm_tfbs_pmcmc_t::options() const {
        return m_options;
}

const data_tfbs_t&
dpm_tfbs_pmcmc_t::data() const {
        return m_data;
}

void
dpm_tfbs_pmcmc_t::m_start_server() {
        vector<save_queue_t<command_t*>* > command_queue;
        for (size_t i = 0; i < m_size; i++) {
                command_queue.push_back(&operator[](i).command_queue());
        }
        if (m_socket_file != "" && m_server == NULL) {
                remove(m_socket_file.c_str());
                m_server = new server_t(m_ios, m_socket_file, command_queue, m_output_queue);
                m_bt     = new boost::thread(boost::bind(&io_service::run, &m_ios));
        }
}

void
dpm_tfbs_pmcmc_t::m_stop_server() {
        if (m_socket_file != "" && m_server != NULL) {
                m_ios.stop();
                m_bt->join();
                delete(m_bt);
                delete(m_server);
                remove(m_socket_file.c_str());
        }
}

void
dpm_tfbs_pmcmc_t::save(const string& filename) const
{
        if (filename == "") {
                cout << *this;
        }
        else {
                ofstream file;
                file.open(filename.c_str());
                file << *this;
                file.close();
        }
}

static
ostream& operator<< (ostream& o, const matrix<double>& m)
{
        for (size_t i = 0; i < m.size(); i++) {
                o << "\t";
                if (m[i].size() == 0)
                        o << "-";
                for (size_t j = 0; j < m[i].size(); j++)
                        o << m[i][j] << " ";
                o << endl;
        }
        return o;
}

ostream& operator<< (ostream& o, const dpm_tfbs_pmcmc_t& pmcmc)
{
        const sampling_history_t& history = pmcmc.sampling_history();

        o << "[Result]" << endl;

        o << "components =" << endl
          << history.components;
        o << "cluster-sizes =" << endl
          << history.cluster_sizes;
        o << "switches =" << endl
          << history.switches;
        o << setprecision(2)
          << fixed;
        o << "likelihood =" << endl
          << history.likelihood;
        o << "posterior =" << endl
          << history.posterior;
        o << "temperature =" << endl
          << history.temperature;
        o << "partitions =" << endl
          << history.partitions;

        return o;
}
