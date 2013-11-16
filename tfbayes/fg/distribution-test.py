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
    y = map(lambda xp: d.density([xp]), x)
    return sum(y)*step

def plot_density(d, x_from, x_to, n=1001):
    x = np.linspace(x_from, x_to, num=n)
    y = map(lambda xp: d.density([xp]), x)
    p = plt.plot(x,y)
    plt.show()

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

# test moments
################################################################################

d = dirac_distribution_t([2])

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
plt.savefig("distribution-test.png")

# distribution on a product space
################################################################################

d = pnormal_distribution_t(2, 1, 2)

delta = 0.05
x = np.arange(-1.0, 3.0, delta)
y = np.arange(-1.0, 3.0, delta)

z = []
for xp in x:
    for yp in y:
        z.append(d.density([xp, yp]))

z = np.array(z)
z = z.reshape(len(x), len(y))

cmap = cm.get_cmap(name='YlOrBr', lut=None)
p = plt.imshow(z, cmap=cmap, origin='lower',
               vmin=z.min(), vmax=z.max(),
               extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plt.show()
