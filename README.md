## Documentation

Please read our [paper](http://arxiv.org/abs/1305.3692) on inference of phylogenetic trees.

## Configuration of local installations

It is necessary to export some environment variables for a local installation. Here is an example:

	export PATH=$HOME/.usr/bin:$PATH
	export CPATH=$HOME/.usr/include
	export LD_LIBRARY_PATH=$HOME/.usr/lib:$LD_LIBRARY_PATH
	export LIBRARY_PATH=$HOME/.usr/lib:$LIBRARY_PATH
	export MANPATH=$HOME/.usr/share/man:$MANPATH

The definitions can for instance be placed in the local *.profile* or *.bash_profile*.

## Requirements

The following libraries are required for tfbayes:

	boost_system
	boost_thread
	boost_regex
	glpk
	gsl
	pthread

Requirements for parsing phylogenetic trees in newick format:

	bison
	flex

Some of the scripts are written in *python* and require:

	biopython
	numpy
	matplotlib
	weblogolib

## Installation

First create all autoconf and automake files with

	autoreconf

For local installations to *$HOME/.usr* use

	./configure --prefix=$HOME/.usr

and otherwise simply

	./configure

The preferred compiler is *clang*, to use it type

	CC=clang CXX=clang++ ./configure

Now the source can be compiled and installed with

	make
	make install

### Known Errors

The macro *__extern_always_inline* may not always be defined and causes an error

	error: unknown type name '__extern_always_inline'

In this case it es necessary to define *CXXFLAGS='-D__extern_always_inline=inline'*.

## Example: Phylogenetic tree inference

The *data* directory contains some data sets for phylogenetic tree inference. We consider the *a subsequence of the MT-RNR2 alignment* for this example. The following command runs *10* Markov chains in parallel, each generates *10000* tree samples:

	tfbayes-treespace-optimize -b 0 -m 10000 -j 10 -p test.posterior.dat metropolis-hastings data/trees/ucsc-hg19-multiz46.nh data/alignments/ucsc-hg19-multiz46-U25123.fa | gzip -f > test.nh.gz

The file *test.posterior.dat* contains the (unnormalized) posterior values for the samples. It is formatted such that it can be easily read and plotted with *GNU-R*. To compute the posterior expectation use

	zcat test.nh.gz > tfbayes-treespace-mean -n 100 -r -f -d 8000 -v mean > test.mean.nh

which drops the first *8000* samples and performs *100* iterations of the algorithm. Similarly, use

	zcat test.nh.gz > tfbayes-treespace-mean -n 100 -r -d 8000 -v median > test.median.nh

to compute the median. To compare the result with the majority rule consensus tree use

	zcat test.nh.gz > tfbayes-treespace-mean -d 8000 -v majority-consensus

## Alignment gaps

The library supports two ways of handling alignment gaps. Which one is used is coded in the alignment data:

+ 'N': The gap is considered as missing data, which means that a nucleotide should be present at this position, but we simply do not know which one (wildcard). It is equivalent to treating the species as if it was not present in the data set, i.e. the species is removed from the phylogenetic tree. An 'N' is commonly used by repeat masking software.
+ '-': The gap is interpreted as an additional character in the alphabet (i.e. a fifth nucleotide). Note that if this is not used in the alignment, the prior counts for this character should be set to zero.

## Newick format

TFBayes uses the following grammar to parse trees in newick format:

	tree_list -> tree_list tree ";"
	tree_list -> tree ";"
	tree      -> "(" node_list "," outgroup ")"
	node_list -> node_list "," node
	node_list -> node
	node      -> "(" node_list "):" distance
	node      -> name ":" distance
	outgroup  -> name ":" distance
	name      -> [a-zA-Z_][a-zA-Z0-9_]*
	distance  -> -?{[0-9]}+("."{[0-9]}*)?

The grammar shows that trees are not required to have binary branching points. However, thee root is expected to have at least three nodes attached to it, i.e. trees are required to have the structure of unrooted trees. The last node attached to the root is required to be a leaf, similar to the convention of MrBayes. Internal edges do not have labels in TFBayes. A valid tree is for instance

	((speTri1:0.322352,cavPor3:0.294901):0.117009,(dipOrd1:0.396332,mm9:0.243578):0.135282,ochPri2:0.340420);

or

	(speTri1:0.322352,cavPor3:0.294901,(dipOrd1:0.396332,mm9:0.243578):0.135282,ochPri2:0.340420);
