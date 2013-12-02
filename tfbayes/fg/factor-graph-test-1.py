#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# example factor graph
################################################################################

def construct_fg(data):
    data = map(lambda x: [x], data)
    fg  = factor_graph_t()
    fg += normal_fnode_t("normal1")
    fg += normal_vnode_t("x")
    fg += normal_fnode_t("normal2", 0, 0.01)
    fg += normal_vnode_t("mu")
    fg += gamma_fnode_t ("gamma", 1, 2)
    fg += gamma_vnode_t ("tau")
    fg.link("normal1:output",    "x")
    fg.link("normal2:output",    "mu")
    fg.link("gamma:output",      "tau")
    fg.link("normal1:mean",      "mu")
    fg.link("normal1:precision", "tau")
    fg.variable_node("x").condition(data)
    return fg

# utility
################################################################################

def plot_density(ax, d, x_limits, xlab="", ylab="", n=1001):
    x = np.linspace(x_limits[0], x_limits[1], num=n)
    y = map(lambda xp: d([xp]), x)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax.plot(x,y)

def plot_fg(fg, data, bound):
    # obtain estimates
    e_mu  = fg["mu" ].moments(0)
    e_tau = fg["tau"].moments(1)
    d     = normal_distribution_t(e_mu, e_tau)
    d1    = fg["mu" ]
    d2    = fg["tau"]
    # plot result
    plt.clf()
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((3,2), (1,0))
    ax3 = plt.subplot2grid((3,2), (1,1))
    ax4 = plt.subplot2grid((3,2), (2,0), colspan=2)
    plot_density(ax1, d,  [ 0.5,  1.5], xlab="x", ylab="density")
    plot_density(ax2, d1, [ 0.9,  1.1], xlab=r'$\mu$', ylab="density")
    plot_density(ax3, d2, [50, 90], xlab=r'$\tau$')
    ax1.hist(data, 50, normed=1)
    ax4.plot(bound)
    ax4.set_xlabel("iteration")
    ax4.set_ylabel("bound")
    #plt.savefig("factor-graph-test-1.png")
    plt.show()

# test
################################################################################

# generate some data
mu    = 1
sigma = 0.1
data  = np.random.normal(mu, sigma, 1000)

# construct and execute the factor graph
fg = construct_fg(data)
bound = fg()

plot_fg(fg, data, bound)
