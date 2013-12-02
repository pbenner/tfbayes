#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# example factor graph
################################################################################

def construct_fg(data):
    data = map(lambda x: [x], data)
    fg  = factor_graph_t()
    fg += categorical_fnode_t("categorical", 3)
    fg += categorical_vnode_t("x", 3)
    fg += dirichlet_fnode_t  ("dirichlet", [2,1,2])
    fg += dirichlet_vnode_t  ("theta", 3)
    fg.link("categorical:output", "x")
    fg.link("categorical:theta",  "theta")
    fg.link("dirichlet:output",   "theta")
    fg.variable_node("x").condition(data)
    return fg

# test
################################################################################

data = [0,0,0,1,1,2,1,1,2,0]

# construct and execute the factor graph
fg = construct_fg(data)
bound = fg()
