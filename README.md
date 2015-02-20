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
The code only works on FASTA files (not FASTQ or any other format).
Here's an example using 10 clusters:
```
julia ARK.jl -i /path/to/FASTA.fa -o /path/to/Output.tsv -n 10
```
Another example using deterministic clustering, 5 clusters, and the SEK training database:
```
julia ARK.jl -i /path/to/FASTA.fa -o /path/to/Output.tsv -n 10 -c Deterministic -t SEK
```
Other options are available, see `julia ARK.jl -h`.

The output format is consistent with the [CAMI challenge](http://www.cami-challenge.org/) and is similar to the output produced by [MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan).

## Further Notes ##
If your installation of dna_utils results in the executable being located in a non-standard location, specify this location using the option ` -k /path/to/./kmer_counts_per_sequence `

It is very important that your installation of BLAS matches the architecture of your hardware (if not, significant increases in computation time might be observed). We recommend using OpenBLAS!