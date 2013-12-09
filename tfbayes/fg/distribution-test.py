#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# utility
################################################################################

def integrate(d, x_from, x_to, n=1001):
    step = (float(x_to) - float(x_from))/(n-1.0)
    x = np.linspace(x_from, x_to, num=n)
    y = map(lambda xp: d([xp]), x)
    return sum(y)*step

def plot_density(d, x_from, x_to, n=1001):
    x = np.linspace(x_from, x_to, num=n)
    y = map(lambda xp: d([xp]), x)
    p = plt.plot(x,y)

def plot_density_2d(d, xlim, ylim, xlab="x", ylab="y", n=101):
    plt.clf()
    ax = plt.subplot(1,1,1)
    x = np.linspace(xlim[0], xlim[1], num=n)
    y = np.linspace(ylim[0], ylim[1], num=n)
    z = []
    for xi in x:
        z.append(map(lambda yp: d([xi, yp]), y))
    p = ax.imshow(z, origin='lower', cmap='Reds', extent=[xlim[0],xlim[1],ylim[0],ylim[1]])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

def plot_density_list(d_list, x_from, x_to, labels=None, xlab="x", ylab="density", n=1001):
    plt.clf()
    p_list = []
    ax = plt.subplot(1,1,1)
    for i, d in enumerate(d_list):
        x = np.linspace(x_from, x_to, num=n)
        y = map(lambda xp: d([xp]), x)
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

d = normal_distribution_t(1,2)

d.moments(0)
d.moments(1)

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
n0.renormalize()

plot_density_list([n0,n1,n2], -5, 5, ["n0", "n1", "n2"])
plt.savefig("distribution-test-1.png")

# bivariate normal distribution
################################################################################

d = binormal_distribution_t(2, 2, 1, 0.9, 0.8)
d([2,2])
plot_density_2d(d, [0,4], [0,4])
plt.savefig("distribution-test-2.png")
