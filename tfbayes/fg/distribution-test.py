#! /usr/bin/env python

from tfbayes.fg import *

import numpy as np

# test if density is normalized
################################################################################

def integrate(f, x_from, x_to, n=1001):
    step = (float(x_to) - float(x_from))/(n-1.0)
    x = np.linspace(x_from, x_to, num=n)
    y = map(f, x)
    return sum(y)*step

integrate(normal_distribution_t(1,2), -5, 5)

# plot density
################################################################################

import matplotlib.pyplot as plt

def plot_density(d, x_from, x_to, n=1001):
    x = np.linspace(x_from, x_to, num=n)
    y = map(d.density, x)
    p = plt.plot(x,y)
    plt.show()

plot_density(normal_distribution_t(1,2), -5, 5)
