#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from tfbayes.fg import *

# utility
################################################################################

def integrate(f, x_from, x_to, n=1001):
    step = (float(x_to) - float(x_from))/(n-1.0)
    x = np.linspace(x_from, x_to, num=n)
    y = map(f, x)
    return sum(y)*step

def plot_density(d, x_from, x_to, n=1001):
    x = np.linspace(x_from, x_to, num=n)
    y = map(d.density, x)
    p = plt.plot(x,y)
    plt.show()

# test if density is normalized
################################################################################

integrate(normal_distribution_t(1,2), -5, 5)

# plot density
################################################################################

plot_density(normal_distribution_t(1,2), -5, 5)

# multiplication test
################################################################################

n1 = normal_distribution_t(1,2)
n2 = normal_distribution_t(2,3)

plot_density(n1, -5, 5)

n1 *= n2

plot_density(n1, -5, 5)

n1.renormalize()

plot_density(n1, -5, 5)
