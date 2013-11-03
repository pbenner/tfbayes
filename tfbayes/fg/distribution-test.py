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
