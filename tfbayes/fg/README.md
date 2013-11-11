
## Factor Graph Library

### Example: Normal data

Construct a factor graph with normal distributed observations, a normal prior for the mean and a gamma prior for the precision parameter.

	def construct_fg1(data):
	    d  = len(data)
	    f1 = pnormal_fnode_t(d, 0, 0, "f1")
	    v1 = data_vnode_t(data, "v1")
	    f2 = normal_fnode_t(0, 0.01, "f2")
	    v2 = normal_vnode_t("v2")
	    f3 = gamma_fnode_t(1, 2, "f3")
	    v3 = gamma_vnode_t("v3")
	    f1.link("output", v1)
	    f2.link("output", v2)
	    f3.link("output", v3)
	    f1.link("mean", v2)
	    f1.link("precision", v3)
	    return factor_graph_t([f1,f2,f3],[v1,v2,v3])


Generate some data and run the factor graph

	mu    = 1
	sigma = 0.1
	data  = np.random.normal(mu, sigma, 1000)

	fg = construct_fg1(data)
	fg()


![alt tag](factor-graph-test-1.png)

### Example: Distributions

Create two normal distributions with different mean and variance. The two distributions are multiplied and the resulting distribution is shown after renormalization.

	n0  = normal_distribution_t()
	n1  = normal_distribution_t(1,2)
	n2  = normal_distribution_t(2,3)
	n0 *= n1
	n0 *= n2
	n0.renormalize()

![alt tag](distribution-test.png)
