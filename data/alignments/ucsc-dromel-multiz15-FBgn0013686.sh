#! /bin/bash

ROOT=ucsc-dromel-multiz15-FBgn0013686

# sample trees
tfbayes-treespace-optimize -b 0 -m 10000 -j 10 -p ${ROOT}.posterior.dat metropolis-hastings ${ROOT}.init.nh ${ROOT}.fa > ${ROOT}.nh

# gzip
gzip ${ROOT}.nh

# compute mean
zcat ${ROOT}.nh.gz | tfbayes-treespace-mean -n 100 -r -f -d 8000 -v mean > ${ROOT}.mean.nh

# compute median
zcat ${ROOT}.nh.gz | tfbayes-treespace-mean -n 100 -r -f -d 8000 -v median > ${ROOT}.median.nh
