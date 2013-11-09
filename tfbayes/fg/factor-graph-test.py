#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from tfbayes.fg import *

# some example factor graphs
################################################################################

def construct_fg1():
    f1 = normal_fnode_t(1, 2, "f1")
    v1 = data_vnode_t(1.42, "v1")
    f2 = normal_fnode_t(2, 2, "f2")
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

def plot_density_list(d_list, x_from, x_to, labels=None, xlab="x", ylab="density", n=1001):
    plt.clf()
    p_list = []
    ax = plt.subplot(1,1,1)
    for i, d in enumerate(d_list):
        x = np.linspace(x_from, x_to, num=n)
        y = map(d.density, x)
        if labels:
            p_list.append(ax.plot(x,y,label=labels[i]))
        else:
            p_list.append(ax.plot(x,y))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if labels:
        ax.legend()

# test
################################################################################

fg = construct_fg1()

fg()

l = ["v2", "v3"]
d = map(lambda name: exponential_family_i(fg[name]), l)

plot_density_list(d, 0, 5, labels=["normal","gamma"])
plt.show()

fg["v2"].moment(1)
fg["v3"].moment(1)
