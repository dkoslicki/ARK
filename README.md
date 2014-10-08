#Warning
*This repository is still in development, please check back soon*

#ARKQuikr

ARKQuikr is an AggRegation of Kmers modification of the QUadratic, Iterative, K-mer based Reconstruction technique  (Quikr) that utilizes sparsity promoting ideas from the field of compressed sensing to reconstruct the composition of a bacterial community (when the input data is a FASTA file of 16S rRNA reads). This extremely fast method comes with a several databases that can be custom trained. Typically reconstruction is accurate down to the genus level.

#What does this repository contain?

This repository is a Julia implementation of the ARKQuikr algorithm. For a Matlab version of this code, see [this website] (http://www.kth.se/en/ees/omskolan/organisation/avdelningar/commth/research/software).


## Requirements ##
+ Mac OS X 10.6.8 or GNU/Linux
+ 4Gb of RAM minimum. Absolutely necessary.
+ gcc that supports OpenMP
+ [dna\_utils](http://github.com/EESI/dna-utils/) must be installed

### Mac Requirements ###
+ Mac OS X 10.6.8 (what we have tested)
+ GCC 4.7 or newer. (gcc 4.2 did not work, and is the default installation)
+ OpenMP libraries (libgomp, usually comes with gcc)

### Linux Requirements ###
+ GCC 4.7 or newer
+ OpenMP libraries (libgomp, usually comes with gcc)

## Installation ##
After cloning and installing the [dna\_utils](http://github.com/EESI/dna-utils/) repository, just clone this repository. As the code contained herein are Julia scripts, no compilation is necessary.


## Usage ##