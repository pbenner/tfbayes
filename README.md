## Documentation

Please read our [paper](http://arxiv.org/abs/1305.3692) on inference of phylogenetic trees.

## Configuration of local installations

It is necessary to export some environment variables for a local installation. Here is an example:

	export PATH=$HOME/.usr/bin:$PATH
	export CPATH=$HOME/.usr/include
	export LD_LIBRARY_PATH=$HOME/.usr/lib:$LD_LIBRARY_PATH
	export LIBRARY_PATH=$HOME/.usr/lib:$LIBRARY_PATH
	export MANPATH=$HOME/.usr/share/man:$MANPATH
	export PYTHONPATH=$HOME/.usr/lib/python2.7/site-packages/:$PYTHONPATH

where of course the system's python version has to be used. The definitions can for instance be placed in the local *.profile* or *.bash_profile*.

## Requirements

The following libraries are required for tfbayes:

	boost (>= 1.47)
	boost_python
	boost_system
	boost_thread
	boost_regex
	glpk
	gsl
	pthread

Requirements for parsing phylogenetic trees in newick format:

	bison (>= 2.7)
	flex

Some of the scripts are written in *python* and require:

	biopython
	numpy
	matplotlib

## Installation

First create all autoconf and automake files with

	autoreconf

For local installations to *$HOME/.usr* use

	./configure --prefix=$HOME/.usr

and otherwise simply

	./configure

The preferred compiler is *clang*, to use it type

	CXX=clang++ ./configure

Now the source can be compiled and installed with

	make
	make install

### Link time optimization (LTO)

LTO can significantly improve the performance of TFBayes. It is recommended to use it with *clang* and *clang++*. It is disabled by default, to switch it on use

	CXX=clang++ ./configure --enable-lto

### Known errors

 * **error: unknown type name '__extern_always_inline'**: The macro *__extern_always_inline* may not always be defined and in this case causes an error. If this happens it es necessary to declare *CXXFLAGS='-D__extern_always_inline=inline'*.
 * **no archive symbol table (run ranlib)**: Your linker is not using the LLVMgold plugin. Either you are not using the gold linker or the plugin is not found.

## Example: Phylogenetic tree inference

The *data* directory contains some data sets for phylogenetic tree inference. We consider the *a subsequence of the MT-RNR2 alignment* for this example. The following command runs *10* Markov chains in parallel, each generates *10000* tree samples:

	tfbayes-treespace-optimize --steps=10000 --chains=10 --save-posterior=test.posterior.dat metropolis-hastings data/trees/ucsc-hg19-multiz46.nh data/alignments/ucsc-hg19-multiz46-U25123.fa | gzip -f > test.nh.gz

The file *test.posterior.dat* contains the (unnormalized) posterior values for the samples. It is formatted such that it can be easily read and plotted with *GNU-R*. To compute the posterior expectation use

	zcat test.nh.gz | tfbayes-treespace-mean -n 100 -r -f -d 8000 -v mean > test.mean.nh

which drops the first *8000* samples and performs *100* iterations of the algorithm. Similarly, use

	zcat test.nh.gz | tfbayes-treespace-mean -n 100 -r -d 8000 -v median > test.median.nh

to compute the median. To compare the result with the majority rule consensus tree use

	zcat test.nh.gz | tfbayes-treespace-mean -d 8000 -v majority-consensus

With *tfbayes-treespace-histogram* several summarizing statistics can be computed. For instance, to obtain a histogram of tree topologies use

	zcat test.nh.gz | tfbayes-treespace-histogram -d 8000 topology > test.topology.dat

which can be visuablized with R

	> attach(read.table("test.topology.dat", header=T))
	> hist(topology)

The topologies are sorted according to their frequencies. With

	zcat test.nh.gz | tfbayes-treespace-histogram -d 8000 edges > test.edges.dat

a table of edge lengths is printed, which can be visualized with

	hist.edges <- function(t, s1, s2, from=-0.2, to=0.2, n=50, main="", ...)
	{
	  x <- c(-t[[s1]], t[[s2]])
	  x <- x[x > from & x < to]
	  hist(x, breaks=seq(from=from,to=to, length.out=n), freq=F,
	       ylab="Density estimate", main=main, ...)
	  lines(density(x,  na.rm=T, adjust=2))
	}
	t <- read.table("test.edges.dat", header=T)
	hist.edges(t, "s14", "s15")

The histogram shows edge lengths of split *s14* as negative values and lengths of split *s15* as positive values. Split identifiers are declared in the header of *test.edges.dat*. After learning a tree it can be used to generate alignments. By comparing generated alignments to some real data one may assess the goodness of the learned tree. For instance, use

	tfbayes-generate-alignment -a 0.2:0.2:0.2:0.2:0.2 simple test.mean.nh

to generate a conserved region

	ochPri2: TCAATACGAG-C-CCACA--GCG--GC-T--GGCTTGCAA-CAA-CA-ATACC-AGTAAGCCGA-TCT-AGGT-G-GACATTCCAGCT-GTCACT-ATAG
	calJac1: TCAAT-CGTGA--CGAC---GCGT-GCTTAGGGCTTACAA-TTA-CTGA-GCC-AAAAGGCCGA-TCACAGCA-GTCACATTC-AGCTCGTAAGA-GTCG
	papHam1: TCAAT-CGAGA--CGACG-CGCGT-GCTT-GGGCTTACAA-TTA-CT-A--CA-CGAAGGCCGA-TCACAGCA-G-GCCACTCC-GCTCGTCACT-ATCG
	panTro2: TCAAT-CGAGA--CGACAG-GCGT-GATT--GGCTTTCAA-CTA-CT-A-ACACCGAAGGCCGA-TCA-ATCA-G-GACATTCCAGCT-GTCACT-ATAG
	   hg19: TCAAT-CGAGA--C-ACA--GCGT-GCTT--GGCTTTCAA-CTA-CT-A-GCACCGAAGGCCGA-TCA-ATCA-G-CCCATTCCAGCT-GTCACT-ATAG
	gorGor1: TCAATCCGAGA--CGACAG-GCGT-GCTT-GGGCTTTCAA-CTA-CT-A-ACACCGAAGGCCGA-TCA-ATCA-G-GCCATTCCAGCTCGTCACT-ATAG
	ponAbe2: TCAAT-CGAGA--CGACA--GCG--GCTT--GGCTTTCAA-CAA-CT-A-ACACCGAAGGCCGA-TCACAGCA-GGGCCATTCCAGCT-GTCACT-ATAG
	    mm9: TCAATACGAGAC-CGAAA--GCGTTGCTT--GGCTATC-A-CAA-CT-G--CACTGAAGGCCGA-TCA-AGCTCG-CACATTCCAGCT-GTCACT-ATAG
	    rn4: TCAATCCGAGAC-CGACA--GCGT-GC-T--GGCTATC-A-CAA-CT-GGACACCGAAGGCCGA-TCACAGCTCG-CACATTC-AGCT-GTCACT-ATTG
	dipOrd1: ACAAT-CGAGAC-CGAAA--GCGT-GC-T--GGCTTTTAG-CAA-CT-GTGCCCCGAACGCCTA-TCA-AGCT-G-GACATTCC-GCT-GT-ACT-ATAG
	cavPor3: ACAAT-CGAG-C-CGACA--G-GT-GC-T--TGCTTTCAA--TA-CTGG-ACACCGAAGGCCGAGTCACCCCT-G-GACAGTCGAGCTCGTCACT-GTAG
	speTri1: TCAAT-CGAGAC-CGACA--GCCT-GCTT--TGCTTACAA-CTA-CT-G-ACACCGAAGGACGA-TCA-AGCT-GGGAAATTCCAGCT-GTCACT-ATAG
	oryCun2: TCAATACGAGAC-CGACA--G-GT-GATT--GGCTTGCAA-CAA-CA-G-ACACAGAACGCCTA-TCT-AGCT-GGGACATTCCAGCT-GTCAC--ATAG

On the other hand, the command

	tfbayes-generate-alignment -a 10:10:10:10:20 simple test.mean.nh

generates a less conserved region with plenty of gaps, i.e.

	ochPri2: --T-G-AG--GC-T-T-CC-CACGAAAGA-CCT-GTA-T-CTCGGTGG-GA--CTGGTGA-AA---C--CCAG-GT--G-CT-G-AA--CAGGATAC--C
	calJac1: C-AACC-CT----T--G-T---CC--AGC-GCC-TT-G--AGAGGG----G--C-GAAGA-TGT-GA-GAGAG-GTT-G-GGCGTAA-TTCAG-CA-G--
	papHam1: C-AACAATT-C--G-GGAT-A-CCA-TGC--C-CTTAG-GA--GGGGG-G--AC-GTAGA-GAT--CCGAGAG-GTT-T-GGTGT-A-TCAA--CA---T
	panTro2: T-AACAGGT-GC-G-GGA--AACCA-TGC-----GT-G-GAT-GGGGG-GAAAC-GGAGA-AAG--CCGCGACGGTT-T-GGAGT-A-TCAT--CA---T
	   hg19: T-AACAGGT-G--G-GGA--AACCA-TGC-----GT-G-GAT-GGGCG-GAAAC-GGAGA-AAG--CC-CGACGGTT-T-GGAGT-A-TCAG--CA---T
	gorGor1: T-CACAGGT-G--G-GGA--AACCA-T-C-----GT-G-GAT-GGGGG-GAAAC-GGAGA-A-G--CCGCGACGGTA-T-GG-GT-A-TCAA--CA---T
	ponAbe2: T--ACAGGT-G--G-GGA--CACCA-TGC--C--TT-G-GAT-GGCAG-GAAAC-GGAGA-AA---CCGCGAC-GTT-T-GGTGT-A-TCAA--CA----
	    mm9: ATAAG-G---G--TGT-TCG-AACAAAGC--C---T-GT-CTCTGCAG-G--ACTGC-GA-AA--GATGCCAG-GT--G-C-T---TC-CAGG-TA--TT
	    rn4: -TAAG-G-T-G--T-TGCCCCAGCAAGGC--C-----GT-CTCAGGAG-GC-ACT-C-GA-AA--G-TGCCAG-GT--G-C-T---G--CAGG-TA--TT
	dipOrd1: A-AA--GGCGG-C--TGCCGATGCTAAGC--C--CT-GTTC-CAA-GG-GA-TCCGGAGA-AA--G-CG-CAG-GT----C-TG-GA--CAG-CTAC-TT
	cavPor3: --AGG-GGT-GC-G--GCC-CA-CTAGCCA-C---T-G--CTCAGGGG----ACTGTAGAGAA--G-CTC-AGGGG--G-C-AGA--A-CAGG-TTC--T
	speTri1: ---GGCG--------T-CG-CAGATAAGC--CC-GT-GT-CTCGG-GG-GG-ACTG-AGACGAA---TACCCC-GT--G-C-TG--AA-CATG-TTC-TT
	oryCun2: ---AGGAGA-G--T-T-CC-CACCGACGC-GCG--CTGT-CTCGGCGA-GA-ACTGGTGA-GA---C-GACAG-GT--G-CT-A-AA-T-AGGATAC-TC


## Example: ChIP-Seq data analysis

Sequences from a ChIP-Seq experiment must be available in *maf* or *mfa* format. In a first step, the training data is preprocessed. For each ChIP-Seq peak in our target species (e.g. DroMel) we are given the nucleotide sequence around this location as a multiple sequence alignment. The purpose of the analysis is to find the motif for our target species and we regard the sequences of all other species as additional information. Therefore, we first remove all columns in the alignment where the target species has a gap ('-'):

	tfbayes-preprocess-alignment -v -s DroMel -m 50 training-set.orig.maf > training-set.filtered.maf

In addition, the command masks all sites in a sequence as missing data ('N') if more than 50 consecutive gaps appear. The filtered data is then used to compute the phylogenetic approximation:

	tfbayes-approximate -v $(PHYLOTREE) < training-set.filtered.maf > training-set.approximation.fa

The sampler requires the alignment data in *fasta* format, we convert the *maf* file with

	tfbayes-maf-to-fasta training-set.filtered.maf training-set.filtered.fa

Before running the sampler, we need to specify a configuration file (*training-set.cfg*):

	[TFBS-Sampler]
	alignment-file    = training-set.filtered.fa
	phylogenetic-file = training-set.approximation.fa
	save = training-set.result
	socket-file = training-set.srv
	process-prior = pitman-yor process
	samples = 1000:100
	population-size = 4
	alpha = 10
	discount = 0.0
	lambda = 0.00000000000001
	initial-temperature = 5
	tfbs-length = 10
	background-model = independence-dirichlet
	background-alpha =
	                 10.0
	                 10.0
	                 10.0
	                 10.0
	                 10.0
	baseline-priors = baseline-default
	baseline-default =
	                 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200
	                 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200
	                 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200
	                 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200 0.200
	                 0.100 0.100 0.100 0.100 0.100 0.100 0.100 0.100 0.100 0.100

The sampler is executed with

	tfbayes-sampler training-set.cfg

which first generates a sequence of 100 burnin samples with temperature greater one and afterwards starts the actual MCMC simulation. The sampler runs 4 Markov chains in parallel, each generates a set of 1000 samples. Once the sampling process has finished, we may plot the posterior probabilities, number of clusters and temperature with

	tfbayes-plot training-set.cfg

A point estimate (i.e. map, mean or median) is computed with the *tfbayes-estimate* command. The computation of the mean and median might take a while and it might be reasonable to only take a subset of the posterior samples, i.e.

	tfbayes-estimate -v --take=1000 mean training-set.cfg

A point estimate can be converted to a logo with

	tfbayes-partition -v -j mean training-set.cfg

which generates a *training-set.pdf* that contains a motif for each cluster (requires *pdftk*).

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
