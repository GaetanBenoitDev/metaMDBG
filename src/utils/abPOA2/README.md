# abPOA: adaptive banded Partial Order Alignment
[![Latest Release](https://img.shields.io/github/release/yangao07/abPOA.svg?label=Release)](https://github.com/yangao07/abPOA/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/abPOA/total.svg?label=Download)](https://github.com/yangao07/abPOA/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abpoa.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/abpoa)
[![PyPI](https://img.shields.io/pypi/dm/pyabpoa.svg?label=pip%20install)](https://pypi.python.org/pypi/pyabpoa)
[![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963)
[![GitHub Issues](https://img.shields.io/github/issues/yangao07/abPOA.svg?label=Issues)](https://github.com/yangao07/abPOA/issues)
[![Build Status](https://img.shields.io/travis/yangao07/abPOA/master.svg?label=Master)](https://travis-ci.org/yangao07/abPOA)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/abPOA/blob/master/LICENSE)
<!-- [![PyPI](https://img.shields.io/pypi/v/pyabpoa.svg?style=flat)](https://pypi.python.org/pypi/pyabpoa) -->
## Updates (v1.4.0)

- Allow to generate multiple consensus sequences from graph alignment (-d)
- Output sequence name by default in MSA format, remove -A/--msa-header

## Getting started
Download the [latest release](https://github.com/yangao07/abPOA/releases):
```
wget https://github.com/yangao07/abPOA/releases/download/v1.4.0/abPOA-v1.4.0.tar.gz
tar -zxvf abPOA-v1.4.0.tar.gz && cd abPOA-v1.4.0
```
Make from source and run with test data:
```
make; ./bin/abpoa ./test_data/seq.fa > cons.fa
```
Or, install via conda and run with test data:
```
conda install -c bioconda abpoa
abpoa ./test_data/seq.fa > cons.fa
```
## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing abPOA via conda](#conda)
  - [Building abPOA from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [To generate one consensus sequence](#gen_1cons)
  - [To generate multiple consensus sequences](#gen_mcons)
  - [To generate row-column multiple sequence alignment](#gen_msa)
  - [To generate graph information in GFA format](#gen_gfa)
  - [To align sequence to an existing graph in GFA/MSA format](#aln_to_gfa)
  - [To generate a plot of the alignment graph](#gen_plot)
- [Input](#input)
- [Output](#output)
  - [Consensus sequence](#cons)
  - [Row-column multiple sequence alignment](#msa)
  - [Full graph information](#gfa)
  - [Plot of alignment graph](#plot)
- [Algorithm description](#description)
  - [Adaptive banding](#banding)
  - [Minimizer-based seeding and partition](#seeding)
  - [Minimizer-based progressive tree](#tree)
  - [Multiple conensus sequences](#mcons)
- [For development](#dev)
- [Evaluation datasets](#eval)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
abPOA is an extended version of [Partial Order Alignment (POA](10.1093/bioinformatics/18.3.452)) 
that performs adaptive banded dynamic programming (DP) with an SIMD implementation. 
abPOA can perform multiple sequence alignment (MSA) on a set of input sequences and 
generate a consensus sequence by applying the [heaviest bundling algorithm](10.1093/bioinformatics/btg109) 
to the final alignment graph.

abPOA can generate high-quality consensus sequences from error-prone long reads and offer 
significant speed improvement over existing tools.

abPOA supports three alignment modes (global, local, extension) and flexible scoring schemes that allow linear, affine and convex gap penalties. 
It right now supports SSE2/SSE4.1/AVX2 vectorization.

For more information, please refer to our [paper](https://dx.doi.org/10.1093/bioinformatics/btaa963) published in Bioinformatics.

## <a name="install"></a>Installation

### <a name="conda"></a>Installing abPOA via conda
On Linux/Unix and Mac OS, abPOA can be installed via
```
conda install -c bioconda abpoa   # install abPOA program
```

### <a name="build"></a>Building abPOA from source files
You can also build abPOA from source files. 
Make sure you have gcc (>=6.4.0) and zlib installed before compiling.
It is recommended to download the [latest release](https://github.com/yangao07/abPOA/releases).
```
wget https://github.com/yangao07/abPOA/releases/download/v1.4.0/abPOA-v1.4.0.tar.gz
tar -zxvf abPOA-v1.4.0.tar.gz
cd abPOA-v1.4.0; make
```
Or, you can use `git clone` command to download the source code.
This gives you the latest version of abPOA, which might be still under development.
```
git clone --recursive https://github.com/yangao07/abPOA.git
cd abPOA; make
```

### <a name="binary"></a>Pre-built binary executable file for Linux/Unix 
If you meet any compiling issue, please try the pre-built binary file:
```
wget https://github.com/yangao07/abPOA/releases/download/v1.4.0/abPOA-v1.4.0_x64-linux.tar.gz
tar -zxvf abPOA-v1.4.0_x64-linux.tar.gz
```

## <a name="usage"></a>General usage
### <a name="gen_1cons"></a>To generate consensus sequence

```
abpoa seq.fa > cons.fa
```

### <a name="gen_mcons"></a>To generate multiple consensus sequences

```
abpoa heter.fa -d2 > 2cons.fa
```

### <a name="gen_msa"></a>To generate row-column multiple sequence alignment in FASTA format

```
abpoa seq.fa -r1 > out.msa
abpoa seq.fa -r2 > out_cons.msa
```

### <a name="gen_gfa"></a>To generate graph information in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format

```
abpoa seq.fa -r3 > out.gfa
```
To include the generated consensus sequence as a path in the GFA file:
```
abpoa seq.fa -r4 > out.gfa
```

### <a name="aln_to_gfa"></a>To align sequence to an existing graph in GFA/MSA format
```
abpoa -i in.gfa seq.fa -r3 > out.gfa
abpoa -i in.msa seq.fa -r1 > out.msa
```
For GFA input file, `S` and `P` lines are required and are used to reconstruct the alignment graph.
For MSA input file, which is generally a FASTA format file, `-` in the sequence indicates the alignment gap.
```
abpoa seq1.fa -r1 > seq1.msa
abpoa -i seq1.msa seq2.fa > cons.fa
```

### <a name="gen_plot"></a>To generate a plot of the alignment graph

```
abpoa seq.fa -g poa.png > cons.fa
```
See [Plot of alignment graph](#plot) for more details about the plot file.

## <a name="input"></a>Input
abPOA works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. The input file is 
expected to contains multiple sequences which will be processed sequentially to perform the iterative 
sequence-to-graph (partial order) alignment.

abPOA can also take a list of filenames as input with option `-l`, where each line is the path to one 
file containing multiple sequences. Each sequence file is then individually aligned by abPOA to generate a
consensus sequence.

## <a name="output"></a>Output
### <a name="cons"></a>Consensus sequence 
By default, abPOA only outputs the consensus sequence generated from the final alignment graph.
It is in FASTA format with the name field set as "Consensus_sequence".
For example:
```
>Consensus_sequence
ACGTGTACACGTTGAC
```

For diploid input sequences, you may want to generate two or more consensus sequences, simply set `-d/--max-num-cons` as a desired value:
```
abpoa heter.fa -d2
```
and this gives you two consensus sequences:
```
>Consensus_sequence_1
CCATTCCCACCATCCTTACCATCAACATCACCATCCCCACCATCCCCAACACCATTCCCACCATCCCTACCATCACCATCACCATCCCCACCAACATCCCCACCACCATCCTCACTACCATCCCCACCACCATTTCCACCATTCCCACCACAGTCACCATCACCCCCACCATCCCCATCATCATCCGCACCATCCCCACCATCCCCACCACCATCTCCATTACCATCCCCACCACCATCTCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCTCGTTACCATCCCCACCACCTTTTCCACCATTCCCACCATCTCCAACACCTCCCCCACCATCATCCCCACCATCCCCACCACCTTCTCCACCATCATTCTCACCATCCCCACCACCATCTCCACCACCATTCTCACCATCTCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCATCATCCCCACCATCC
>Consensus_sequence_2
CCATTCCCACCATCCTTACCATCAACATCACCATCCCCACCATCCCCAACACCATTCCCACCATCCCTACCATCACCATCACCATCCCCACCAACATCCCCACCACCATCCTCACTACCATCCCCACCACCATTTCCACCATTCCCACCACAGTCACCATCACCCCCACCATCCCCATCATCATCCGCACCATCCCCACCATCCCCACCACCATCTCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCTCGTTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCTCGTTACCATCCCCACCACCTTTTCCACCATTCCCACCATCATCCCCACCGCCATCCTCGTTACCATCCCCACCACCTTTTCCACCATTCCCACCATCTCCAACACCTCCCCCACCATCATCCCCACCATCCCCACCACCTTCTCCACCATCATTCTCACCATCCCCACCACCATCTCCACCACCATTCTCACCATCTCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCATCATCCCCACCATCC
```
### <a name="msa"></a>Row-column multiple sequence alignment
abPOA can also output the row-column multiple sequence alignment (RC-MSA) of all the aligned sequences in FASTA format.
For example:
```
>1
ACGTGTACA-GTTGAC
>2
A-G-GTACACGTT-AC
>3
A-GTGT-CACGTTGAC
>4
ACGTGTACA--TTGAC
```
The `-` in the sequence stands for alignment gap. 

### <a name="gfa"></a>Full graph information
abPOA can output the final alignment graph in GFA format.
Each segment line (`S` line) represents one node and each link line (`L` line) represents one edge between two nodes.
The original input sequences and the generated consensus sequence are described as paths in `P` lines.

abPOA outputs two graph-related numbers in the header line (`H` line):
`NS` and `NL`, which denote the total number of nodes and edges in the GFA file, respectively.

Please refer to the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for more details of the GFA format.

### <a name="plot"></a>Plot of alignment graph

abPOA can generate a plot of the final partial order alignment graph with the help of `graphviz dot`. 
For example:

![pog](https://github.com/yangao07/abPOA/blob/master/pog.png)

The numbers inside the nodes are the node IDs. The numbers on the edges are the edge weights.
`S` and `E` are the auxiliary start and end nodes that have no sequence bases.

Make sure you have `dot` installed before using abPOA to generate the plot.
For Linux/Unix systems: `sudo apt-get install graphviz`.

## <a name="description"></a>Algorithm description
### <a name="banding"></a>Adaptive banding
To understand how the adaptive banding working, please refer to our [Bioinformatics paper](https://dx.doi.org/10.1093/bioinformatics/btaa963).

### <a name="seeding"></a>Minimizer-based seeding mode
As abPOA always allocates quadratic size of memory, for very long input sequences (>10 kb), memory usage will be a challenge.

To solve this issue, we develop a minimizer-based seeding and partition method to split the sequence and graph with a small window.
The full POA DP matrix can be split into several smaller ones and adaptive banded POA can be performed within each small window separately.

In more detail, abPOA extracts all the minimizers from all the input sequences, then all the minimizer hits between each pair of two sequences can be found.
For each pair of sequences, the minimizer hits are first chained together using relatively stringent criteria to make sure that no big gap exists in the chain.
This usually leads to several separated local chains of minimizer hits.
A second round of chaining is then performed on all the local minimizer chains to generate a global chain going through the entire sequence.
With this global chain, abPOA selects a series of minimizer hits as partition anchors which has at least a distance of 500 bp (by default, -n/--min-poa-win).
Within each partitioned window, abPOA performs banded partial order alignment separately and combines all the alignment results at the end.

### <a name="tree"></a>Minimizer-based progressive tree
Instead of aligning all the sequences in the original order, abPOA can alternatively build a progressive tree to guide the alignment order.
The generation of the progressive tree is also based on minimizers.
For each pair of sequences, abPOA calculates their similarity score which is the Jaccard similarity of the minimizers, i.e. the number of minimizer hits divided by the total number of all minimizers from the two sequences.
With all the similarity scores (minimizer-based Jaccard similarity), abPOA builds the progressive tree in the following way:

1. Pick the first two sequences that have the highest scores. The progressive tree set is initialized as these first two sequences.
2. For each remaining sequence, sum the scores between the remaining sequence and all the sequences from the current progressive tree set. Pick the one with the highest sum score, and push it to the progressive tree set.
3. Repeat step 2, until no sequence remains.

Then, abPOA performs partial order alignment following the order of sequences in this progressive tree set.

### <a name="mcons"></a>Multiple consensus sequences
Since v1.4.0, abPOA supports generating multiple consensus sequences from the final alignment graph (set -d/--max-num-cons as >1).

The general underlying idea is to group input sequences into multiple clusters based on the heterozygous bases in the graph,
Then, one consensus sequence is separately generated for each cluster of input sequences.
The minimum allele frequency for each heterozygous base is 0.25 (by default, -q/--min-freq). 

## <a name="dev"></a>For development
abPOA is not only a stand-alone tool for MSA and consensus calling, it can also work as a programming library. [example.c](example.c) shows how to use the C APIs of abPOA to take a set of sequences as input and perform MSA and consensus calling. Basically, the library file `libabpoa.a` and two header files [abpoa.h](include/abpoa.h) and [simd_instruction.h](include/simd_instruction.h) are needed to make the abPOA library work in your program.

abPOA also provides Python bindings to all the primary C APIs. Refer to [python/README.md](python/README.md) for more details.

## <a name="eval"></a>Evaluation datasets
The evaluation datasets and scripts used in [abPOA paper](https://dx.doi.org/10.1093/bioinformatics/btaa963) can be found in [abPOA-v1.0.5](https://github.com/yangao07/abPOA/releases/tag/v1.0.5).

## <a name="contact"></a>Contact
Yan Gao gaoy286@mail.sysu.edu.cn

Yi Xing xingyi@email.chop.edu

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/abPOA/issues)
