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

### Link time optimization (LTO)

LTO can significantly improve the performance of TFBayes. It is recommended to use it with *clang* and *clang++*. It is disabled by default, to switch it on use

	CC=clang CXX=clang++ ./configure --enable-lto

### Known errors

 * **error: unknown type name '__extern_always_inline'**: The macro *__extern_always_inline* may not always be defined and in this case causes an error. If this happens it es necessary to declare *CXXFLAGS='-D__extern_always_inline=inline'*.
 * **no archive symbol table (run ranlib)**: Your linker is not using the LLVMgold plugin. Either you are not using the gold linker or the plugin is not found.

## Example: Phylogenetic tree inference

The *data* directory contains some data sets for phylogenetic tree inference. We consider the *a subsequence of the MT-RNR2 alignment* for this example. The following command runs *10* Markov chains in parallel, each generates *10000* tree samples:

	tfbayes-treespace-optimize -b 0 -m 10000 -j 10 -p test.posterior.dat metropolis-hastings data/trees/ucsc-hg19-multiz46.nh data/alignments/ucsc-hg19-multiz46-U25123.fa | gzip -f > test.nh.gz

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
