#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from tfbayes.fg import *

# example factor graph
################################################################################

def construct_fg(data):
    data = map(lambda x: [x], data)
    m  =  mixture_fnode_t("m1")
    m  += normal_fnode_t("f1")
    m  += normal_fnode_t("f2")
    fg  = factor_graph_t()
    fg += m
    fg += normal_vnode_t("x")
    fg += normal_vnode_t("vmu1")
    fg += normal_vnode_t("vmu2")
    fg += normal_fnode_t("fmu1", 0, 10)
    fg += normal_fnode_t("fmu2", 0, 10)
    fg += categorical_fnode_t("fz", 2)
    fg += categorical_vnode_t("vz", 2)
    fg += dirichlet_fnode_t("fdir", [1,1])
    fg += dirichlet_vnode_t("vdir", 2)
    fg.link("fz:output", "vz")
    fg.link("fmu1:output", "vmu1")
    fg.link("fmu2:output", "vmu2")
    fg.link("m1:indicator", "vz")
    fg.link("m1:f1:mean", "vmu1")
    fg.link("m1:f2:mean", "vmu2")
    fg.link("m1:output", "x")
    fg.link("fdir:output", "vdir")
    fg.link("fz:theta", "vdir")
    fg.variable_node("x").condition(data)
    return fg

# test
################################################################################

data = [0]

# construct and execute the factor graph
fg = construct_fg(data)
bound = fg()
