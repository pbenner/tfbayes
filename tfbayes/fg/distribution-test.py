#! /usr/bin/env python

from tfbayes.fg import *

import numpy as np

normal = normal_distribution_t(0,2)

# test if density is normalized
################################################################################

import matplotlib.pyplot as plt

n      = 1001
x_from = -5
x_to   =  5
step   = (float(x_to) - float(x_from))/(n-1.0)

x = np.linspace(x_from, x_to, num=n)
y = map(normal.density, x)
sum(y)*step

p = plt.plot(x,y)
plt.show()
