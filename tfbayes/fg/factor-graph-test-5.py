#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# example factor graph
################################################################################

def construct_fg(data):
    categorical_vnode_t("z", 3)
    categorical_vnode_t("z", 3)
    # prepare the mixture node
    m   = mixture_fnode_t("mixture")
    m  += normal_fnode_t("normal1", 0, 100)
    m  += normal_fnode_t("normal2", 0, 100)
    # construct graph inside the plate
    fg  = factor_graph_t()
    fg += m
    fg += normal_vnode_t("x")
    fg += categorical_fnode_t("categorical", 2)
    fg += categorical_vnode_t("z", 2)
    fg.link("mixture:indicator",  "z")
    fg.link("mixture:output",     "x")
    fg.link("categorical:output", "z")
    fg.replicate(len(data)-1)
    # construct the remaining graph
    fg += normal_vnode_t("mu1")
    fg += normal_vnode_t("mu2")
    fg += normal_fnode_t("normal1", 0.0, 0.01)
    fg += normal_fnode_t("normal2", 0.0, 0.01)
    fg += dirichlet_fnode_t("dirichlet", [1,1])
    fg += dirichlet_vnode_t("theta", 2)
    fg.link("dirichlet:output",   "theta")
    fg.link("categorical:theta",  "theta")
    fg.link("normal1:output", "mu1")
    fg.link("normal2:output", "mu2")
    fg.link("mixture:normal1:mean", "mu1")
    fg.link("mixture:normal2:mean", "mu2")
    for i, x in enumerate(data):
        fg.variable_node("x", i).condition([[x]])
    return fg

# utility
################################################################################

def plot_density(ax, d, x_limits, xlab="", ylab="", n=1001):
    x = np.linspace(x_limits[0], x_limits[1], num=n)
    y = map(lambda xp: d([xp]), x)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax.plot(x,y)

def plot_fg(fg, data1, data2, bound):
    # obtain estimates
    e_mu1 = fg["mu1"].moments(0)
    e_mu2 = fg["mu2"].moments(0)
    e_d1  = normal_distribution_t(e_mu1, 100)
    e_d2  = normal_distribution_t(e_mu2, 100)
    if e_mu1 < e_mu2:
        d1 = fg["mu1"]
        d2 = fg["mu2"]
    else:
        d1 = fg["mu2"]
        d2 = fg["mu1"]
    # plot result
    plt.clf()
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((3,2), (1,0))
    ax3 = plt.subplot2grid((3,2), (1,1))
    ax4 = plt.subplot2grid((3,2), (2,0), colspan=2)
    plot_density(ax1, e_d1, [-0.5, 0.5], xlab="x", ylab="density")
    plot_density(ax1, e_d2, [-0.5, 0.5], xlab="x", ylab="density")
    plot_density(ax2, d1, [-0.2, 0.0], xlab=r'$\mu_1$', ylab="density")
    plot_density(ax3, d2, [ 0.0, 0.2], xlab=r'$\mu_2$', ylab="density")
    ax1.hist(data1, 50, normed=1)
    ax1.hist(data2, 50, normed=1)
    ax4.plot(bound)
    ax4.set_xlabel("iteration")
    ax4.set_ylabel("bound")
    plt.tight_layout()
    #plt.savefig("factor-graph-test-5.png")
    plt.show()

# test 1
################################################################################

def test1():
    data = [0.0,1.0]
    fg = construct_fg(data)
    print fg()
    for i in range(len(data)):
        print "z[{}]: {{{}, {}}}".format(i, fg["z",i].moments(0), fg["z",i].moments(1))
    print "mu1: ", fg["mu1"].moments(0)
    print "mu2: ", fg["mu2"].moments(0)

# test 2
################################################################################

def test2():
    # generate some data
    data1 = np.random.normal(-0.1, 0.1, 100)
    data2 = np.random.normal( 0.1, 0.1, 200)
    data  = np.append(data1, data2)
    # construct and execute the factor graph
    fg = construct_fg(data)
    bound = fg()
    plot_fg(fg, data1, data2, bound[3:])

# test
################################################################################

#test1()
test2()
