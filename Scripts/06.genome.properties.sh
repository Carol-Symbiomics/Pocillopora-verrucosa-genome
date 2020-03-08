#!/bin/bash

# Genome properties K-hmer based using jellyfish

#Kmer 25
jellyfish count -t 20 -C -m 25 -s 5G -o Pver_25mer *.fq

jellyfish histo -h 3000000 -t 20 -o Pver_25mer.histo Pver_25mer

#Kmer 31
jellyfish count -t 20 -C -m 31 -s 5G -o Pver_31mer *.fq

jellyfish histo -h 3000000 -t 20 -o Pver_31mer.histo Pver_31mer

# Plots and kmer count was done using an the R script Pver.jellyfish.khmer.spectra.R in the Supplementary_scripts folder
