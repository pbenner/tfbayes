#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# some example factor graphs
################################################################################

# link every data point to a single node
def construct_fg0(data):
    fg  = factor_graph_t()
    # factor graph inside the plate
    fg += normal_fnode_t("f1", 0, 0)
    fg += data_vnode_t  ("v1")
    fg.link("f1", "output",    "v1")
    # replicate this graph n-1 times
    fg.replicate(len(data)-1)
    # construct graph outside the plate
    fg += normal_fnode_t("f2", 0, 0.01)
    fg += normal_vnode_t("v2")
    fg += gamma_fnode_t ("f3", 1, 2)
    fg += gamma_vnode_t ("v3")
    fg.link("f2", "output",    "v2")
    fg.link("f3", "output",    "v3")
    # connect v2 and v3 to all factors f1
    # within the plate
    fg.link("f1", "mean",      "v2")
    fg.link("f1", "precision", "v3")
    # loop over all variable nodes v1
    for i, d in enumerate(data):
        fg.data_vnode("v1", i).condition([d])
    return fg

# use a product normal distribution
def construct_fg1(data):
    fg  = factor_graph_t()
    fg += pnormal_fnode_t("f1", len(data), 0, 0)
    fg += data_vnode_t   ("v1")
    fg += normal_fnode_t ("f2", 0, 0.01)
    fg += normal_vnode_t ("v2")
    fg += gamma_fnode_t  ("f3", 1, 2)
    fg += gamma_vnode_t  ("v3")
    fg.link("f1", "output",    "v1")
    fg.link("f2", "output",    "v2")
    fg.link("f3", "output",    "v3")
    fg.link("f1", "mean",      "v2")
    fg.link("f1", "precision", "v3")
    fg.data_vnode("v1").condition(data)
    return fg

# utility
################################################################################

def plot_density(ax, d, x_limits, xlab="", ylab="", n=1001):
    x = np.linspace(x_limits[0], x_limits[1], num=n)
    y = map(lambda xp: d.density([xp]), x)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax.plot(x,y)

# test 1
################################################################################

fg = construct_fg1([3.0])

fg()
fg["v2"].moment(1)
fg["v3"].moment(1)

fg = construct_fg0([3.0])

fg()
fg["v2"].moment(1)
fg["v3"].moment(1)

# test 2
################################################################################

# generate some data
mu    = 1
sigma = 0.1
data  = np.random.normal(mu, sigma, 1000)

plt.hist(data, 50, normed=1)
plt.show()

# construct and execute the factor graph
fg = construct_fg1(data)
fg()

# obtain estimates
e_mu  = fg["v2"].moment(1)[0]
e_tau = fg["v3"].moment(1)[0]
d     = normal_distribution_t(e_mu, e_tau)
d1    = exponential_family_i(fg["v2"])
d2    = exponential_family_i(fg["v3"])

# plot result
plt.clf()
ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
ax2 = plt.subplot2grid((2,2), (1,0))
ax3 = plt.subplot2grid((2,2), (1,1))
plot_density(ax1, d,  [ 0.5,  1.5], xlab="x", ylab="density")
plot_density(ax2, d1, [ 0.9,  1.1], xlab=r'$\mu$', ylab="density")
plot_density(ax3, d2, [50, 90], xlab=r'$\tau$')
ax1.hist(data, 50, normed=1)
#plt.savefig("factor-graph-test-1.png")
plt.show()
