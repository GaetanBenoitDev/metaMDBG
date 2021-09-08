# Bloocoo 

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Bloocoo/job/tool-bloocoo-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Bloocoo/job/tool-bloocoo-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Bloocoo/job/tool-bloocoo-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Bloocoo/job/tool-bloocoo-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

#What is Bloocoo?
Bloocoo is a k-mer spectrum-based read error corrector, designed to correct large datasets with a very low memory footprint. It uses the disk streaming k-mer counting algorithm contained in the GATB library, and inserts solid k-mers in a bloom-filter. The correction procedure is similar to the Musket multistage approach.  Bloocoo yields similar results while requiring far less memory: as an example, it can correct whole human genome re-sequencing reads at 70 x coverage with less than 4GB of memory.

G. Benoit, D. Lavenier, C. Lemaitre, G. Rizk.  (2015) [Bloocoo, a memory efficient read corrector](https://hal.inria.fr/hal-01092960). Inria-HAL.
								
# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions

    # get a local copy of source code
    git clone --recursive https://github.com/GATB/bloocoo.git
    
    # compile the code an run a simple test on your computer
    cd bloocoo
    sh INSTALL

# User manual	 
								
## Description

Bloocoo is a kmer-spectrum based  read error corrector. In a first pass, all  k-mers are counted, then  k-mers more abundant than a given threshold are kept, i.e. “solid k-mers”.

Correction is then performed by scanning  k-mers of a read. For example, a single isolated error generates a gap of k non solid k-mers making the detection of its exact location easy. Correction is made by trying the three different possible nucleotides at the error site, and checking if corresponding k-mers are in the set of solid k-mers.

When several close errors occurs, the pattern is more complex, errors are corrected via a vote algorithm similar to the one in the Musket software (http://musket.sourceforge.net/).

What makes Bloocoo different is the k-mer counting stage and the way solid k-mers are stored in memory. k-mer counting is conducted via the DSK algorithm included in the GATB library, which requires constant-memory. Solid k-mers are stored in a Bloom filter which is fast and memory-efficient : we use only 11 bits of memory per solid k-mers. Therefore,  correction of a whole human genome sequencing read set needs only 4GB of memory.


## Usage

A typical command line is:

    Bloocoo -file reads.fasta -kmer-size 27  -abundance 4

There is 1 mandatory argument:

    -file : the read file name, can be fasta, fastq, gzipped or not.

Two important arguments:

    -kmer-size : the k-mer size (typically ~31)
    -abundance-min : the minimal abundance threshold defining solid k-mers (typically  between 3 and 6, but depends on the read depth, you can also use 'auto' and it is automatically inferred from the data)

Additional useful options :

    -nb-cores  : number of threads used
    -high-recall  :  correct more errors but can also introduce more mistakes
    -slow : slower modes with more pass, but better correction
    -high-precision :  correct safely, correct less errors but introduce less mistakes
    -ion : (experimental) mode for correcting indels present in  ion torrent reads


## Examples

    ./Bloocoo -file reads.fasta
    -> generates the file reads_corrected.fasta
 

Note : 
In order to use k values larger than 31, recompilation is necessary (for the moment, this will be improved in next versions).

In the sequence of commands given in the INSTALL file, change the command: 

    cmake ..

by 

    cmake -DKSIZE_LIST="64" ..

this will allow to use k<63

For larger k, change the value such that it is a multiple of 32

#Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/