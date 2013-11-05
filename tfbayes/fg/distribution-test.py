#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from tfbayes.fg import *

# utility
################################################################################

def integrate(d, x_from, x_to, n=1001):
    step = (float(x_to) - float(x_from))/(n-1.0)
    x = np.linspace(x_from, x_to, num=n)
    y = map(d.density, x)
    return sum(y)*step

def plot_density(d, x_from, x_to, n=1001):
    x = np.linspace(x_from, x_to, num=n)
    y = map(d.density, x)
    p = plt.plot(x,y)
    plt.show()

def plot_density_list(d_list, x_from, x_to, labels=None, xlab="x", ylab="density", n=1001):
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

# test moments
################################################################################

d = dirac_distribution_t(2)
d.moment(1)
d.moment(2)

# test if density is normalized
################################################################################

integrate(normal_distribution_t(1,2), -5, 5)
integrate(gamma_distribution_t(1,2), 0.0, 5)

# plot density
################################################################################

plot_density(normal_distribution_t(1,2), -5, 5)
plot_density(gamma_distribution_t(1,2), 0, 5)

# multiplication test
################################################################################

n1  = normal_distribution_t(1,2)
n2  = normal_distribution_t(2,3)
n1 *= n2

n1.renormalize()

plot_density(n1, -5, 5)

g1  = gamma_distribution_t(1.2,2.3)
g2  = gamma_distribution_t(2,3)
g1 *= g2

g1.renormalize()

plot_density(g1, 0, 5)

# make some nice graphics
################################################################################

n0  = normal_distribution_t()
n1  = normal_distribution_t(1,2)
n2  = normal_distribution_t(2,3)
n0 *= n1
n0 *= n2

plot_density_list([n0,n1,n2], -5, 5, ["n0", "n1", "n2"])
plt.savefig("distribution-test.png")
