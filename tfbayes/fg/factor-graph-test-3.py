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
    fg += categorical_fnode_t("f1", [0,0,0])
    fg += categorical_data_t ("v1", 3)
    fg += dirichlet_fnode_t  ("f2", [2,1,2])
    fg += dirichlet_vnode_t  ("v2", 3)
    fg.link("f1:output", "v1")
    fg.link("f1:theta",  "v2")
    fg.link("f2:output", "v2")
    fg.variable_node("v1").condition(data)
    return fg

# test
################################################################################

data = [0,0,0,1,1,2,1,1,2,0]

# construct and execute the factor graph
fg = construct_fg(data)
bound = fg()
