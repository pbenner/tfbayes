#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# some example factor graphs
################################################################################

# link every data point to a single node
def construct_fg0(data):
    d  = len(data)
    f2 = normal_fnode_t(0, 0.01, "f2")
    v2 = normal_vnode_t("v2")
    f3 = gamma_fnode_t(1, 2, "f3")
    v3 = gamma_vnode_t("v3")
    f2.link("output", v2)
    f3.link("output", v3)
    fnodes = [f2, f3]
    vnodes = [v2, v3]
    for x in data:
        f1 = normal_fnode_t(0, 0, "f1")
        v1 = data_vnode_t([x], "v1")
        f1.link("output", v1)
        f1.link("mean", v2)
        f1.link("precision", v3)
        fnodes.append(f1)
        vnodes.append(v1)
    return factor_graph_t(fnodes, vnodes)

# use a product normal distribution
def construct_fg1(data):
    d  = len(data)
    f1 = pnormal_fnode_t(d, 0, 0, "f1")
    v1 = data_vnode_t(data, "v1")
    f2 = normal_fnode_t(0, 0.01, "f2")
    v2 = normal_vnode_t("v2")
    f3 = gamma_fnode_t(1, 2, "f3")
    v3 = gamma_vnode_t("v3")
    f1.link("output", v1)
    f2.link("output", v2)
    f3.link("output", v3)
    f1.link("mean", v2)
    f1.link("precision", v3)
    return factor_graph_t([f1,f2,f3],[v1,v2,v3])

# utility
################################################################################

def plot_density(ax, d, x_limits, xlab="", ylab="", n=1001):
    x = np.linspace(x_limits[0], x_limits[1], num=n)
    y = map(lambda xp: d.density([xp]), x)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax.plot(x,y)

def plot_density2d(ax, d1, d2, x_limits, y_limits, xlab="", ylab="", n=101):
    x = np.arange(x_limits[0], x_limits[1], float(x_limits[1]-x_limits[0])/n)
    y = np.arange(y_limits[0], y_limits[1], float(y_limits[1]-y_limits[0])/n)
    print y
    z = []
    for xp in x:
        for yp in y:
            z.append(d1.density([xp])*d2.density([yp]))
    # reshape the array
    z = np.array(z)
    z = z.reshape(len(x), len(y))
    # plot density
    cmap = cm.get_cmap(name='YlOrBr', lut=None)
    p = ax.imshow(z, cmap=cmap, origin='lower',
                  vmin=z.min(), vmax=z.max(),
                  extent=[x.min(), x.max(), y.min(), y.max()],
                  aspect='auto')
    return p

def plot_density_list(d_list, x_from, x_to, labels=None, xlab="x", ylab="density", n=1001):
    plt.clf()
    p_list = []
    ax = plt.subplot(1,1,1)
    for i, d in enumerate(d_list):
        x = np.linspace(x_from, x_to, num=n)
        y = map(lambda xp: d.density([xp]), x)
        if labels:
            p_list.append(ax.plot(x,y,label=labels[i]))
        else:
            p_list.append(ax.plot(x,y))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if labels:
        ax.legend()

def plot_fg(fg, nodes, x_from, x_to, labels=None, xlab="x", ylab="density"):
    densities = map(lambda name: exponential_family_i(fg[name]), nodes)
    plot_density_list(densities, x_from, x_to, labels, xlab=xlab, ylab=ylab)

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
plt.savefig("factor-graph-test-1.png")
plt.show()