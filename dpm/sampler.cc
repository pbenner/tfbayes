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

#include <gsl/gsl_randist.h>

#include <sampler.hh>

using namespace std;

GibbsSampler::GibbsSampler(DPM& dpm, const Data& data)
        : _dpm(dpm), _data(data), _sampling_steps(0),
          _sampling_history(*new sampling_history_t())
{
        // for sampling statistics
        _sampling_history.switches.push_back(vector<double>());
        _sampling_history.likelihood.push_back(vector<double>());
        _sampling_history.components.push_back(vector<size_t>());
        _sampling_history.switches[0].push_back(0);
        _sampling_history.likelihood[0].push_back(_dpm.likelihood());
        _sampling_history.components[0].push_back(0);
}

GibbsSampler::GibbsSampler(const GibbsSampler& sampler)
        : _dpm(*sampler._dpm.clone()), _data(sampler._data), _sampling_steps(0),
          _sampling_history(*new sampling_history_t(sampler._sampling_history))
{
}

GibbsSampler::~GibbsSampler()
{
        delete(&_dpm);
        delete(&_sampling_history);
}

GibbsSampler*
GibbsSampler::clone() const {
        return new GibbsSampler(*this);
}

bool
GibbsSampler::_sample(const element_t& element) {
        const word_t word = _data.get_word(element, _dpm.word_length());
        ////////////////////////////////////////////////////////////////////////
        // check if we can sample this element
        if (!_dpm.valid_for_sampling(element, word)) {
                return false;
        }
        ////////////////////////////////////////////////////////////////////////
        // release the element from its cluster
        cluster_tag_t old_cluster_tag = _dpm.cluster_manager().get_cluster_tag(element);
        _dpm.remove_word(word, old_cluster_tag);
        size_t components = _dpm.mixture_components();
        double weights[components+1];
        cluster_tag_t tags[components+1];
        _dpm.mixture_weights(word, weights, tags);

        ////////////////////////////////////////////////////////////////////////
        // draw a new cluster for the element and assign the element
        // to that cluster
        gsl_ran_discrete_t* gdd  = gsl_ran_discrete_preproc(components+1, weights);
        cluster_tag_t i = gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);
        cluster_tag_t new_cluster_tag = tags[i];

        ////////////////////////////////////////////////////////////////////////
        _dpm.add_word(word, new_cluster_tag);

        return old_cluster_tag != new_cluster_tag;
}

const DPM&
GibbsSampler::model() const {
        return _dpm;
}

const sampling_history_t&
GibbsSampler::sampling_history() const {
        return _sampling_history;
}

const posterior_t&
GibbsSampler::posterior() const {
        return _dpm.posterior();
}

void
GibbsSampler::sample(size_t n, size_t burnin) {
        // burn in sampling
        for (size_t i = 0; i < burnin; i++) {
                printf("Burn in... [%u][Components: %02d]\n", (unsigned int)i+1, (int)_dpm.mixture_components());
                fflush(stdout);
                double sum = 0;
                for (Data::const_iterator_randomized it = _data.begin_randomized();
                     it != _data.end_randomized(); it++) {
                        const bool switched =_sample(**it);
                        if (switched) sum+=1;
                }
                _sampling_history.likelihood[0].push_back(_dpm.likelihood());
                _sampling_history.components[0].push_back(_dpm.mixture_components());
                _sampling_history.switches[0].push_back(sum/(double)_data.size());
        }
        // sample `n' times
        for (size_t i = 0; i < n; i++) {
                // loop through all elements
                printf("Sampling... [%u][Components: %02d]\n", (unsigned int)i+1, (int)_dpm.mixture_components());
                fflush(stdout);
                double sum = 0;
                for (Data::const_iterator_randomized it = _data.begin_randomized();
                     it != _data.end_randomized(); it++) {
                        const bool switched = _sample(**it);
                        if (switched) sum+=1;
                }
                _sampling_history.likelihood[0].push_back(_dpm.likelihood());
                _sampling_history.components[0].push_back(_dpm.mixture_components());
                _sampling_history.switches[0].push_back(sum/(double)_data.size());
                _dpm.update_posterior(_sampling_steps);
                _sampling_steps++;
        }
}
