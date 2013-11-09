
from tfbayes.fg import *

# some example factor graphs
################################################################################

def construct_fg1():
    f1 = normal_fnode_t(1, 2, "f1")
    v1 = data_vnode_t(1.42, "v1")
    f2 = normal_fnode_t(2, 2, "f2")
    v2 = normal_vnode_t("v2")
    f3 = gamma_fnode_t(1, 2, "f3")
    v3 = gamma_vnode_t("v3")
    f1.link("output", v1)
    f2.link("output", v2)
    f3.link("output", v3)
    f1.link("mean", v2)
    f1.link("precision", v3)
    return factor_graph_t([f1,f2,f3],[v1,v2,v3])

# test
################################################################################

fg = construct_fg1()

fg()

fg["v2"].moment(1)
fg["v3"].moment(1)
