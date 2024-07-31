MetaMDBG is a fast and low-memory assembler for long and accurate metagenomics reads (e.g. PacBio HiFi). It is based on the [minimizer de-Brujin graph](https://github.com/ekimb/rust-mdbg) (MDBG), which have been reimplemetend specifically for metagenomics assembly. MetaMDBG combines an efficient multi-k approach in minimizer-space for dealing with uneven species coverages, and a novel abundance-based filtering method for simplifying strain complexity.

Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

## Installation

### Conda

```sh
conda install -c conda-forge -c bioconda metamdbg
```
### Building from source (using conda)

<details><summary>See details</summary>
<p>
Choose an installation directory, then copy-paste the following commands.
	
```sh
# Download metaMDBG repository  
git clone https://github.com/GaetanBenoitDev/metaMDBG.git

# Create metaMDBG conda environment
cd metaMDBG
conda env create -f conda_env.yml
conda activate metamdbg
conda env config vars set CPATH=${CONDA_PREFIX}/include:${CPATH}
conda deactivate

# Activate metaMDBG environment
conda activate metamdbg

# Compile the software
mkdir build
cd build
cmake ..
make -j 3
conda install -c bioconda -c conda-forge metamdbg
```
	
After successful installation, an executable named metaMDBG will appear in ./build/bin.
</p>
</details>

### Building from source

<details><summary>See details</summary>
	
<p>
	
**Prerequisites**
- gcc 9.4+
- cmake 3.10+
- zlib
- openmp
- minimap2 2.24+
- wfmash
- samtools 1.6+ (using htslib)
  
</p>
	
```sh
git clone https://github.com/GaetanBenoitDev/metaMDBG.git
cd metaMDBG
mkdir build
cd build
cmake ..
make -j 3
```

</details>

## Usage

```sh
./metaMDBG asm outputDir reads... {OPTIONS}

	outputDir     # Output dir for contigs and temporary files
	reads...      # Read filename(s) (separated by space)
	-t            # Number of cores [3]
	
Examples:
./metaMDBG asm ./path/to/assemblyDir reads.fastq.gz -t 4                                        #single-sample assembly
./metaMDBG asm ./path/to/assemblyDir reads_A.fastq.gz reads_B.fastq.gz reads_C.fastq.gz -t 4    #co-assembly
```

MetaMDBG will generate polished contigs in outputDir ("contigs.fasta.gz").

## Input data (PacBio and Nanopore)

- MetaMDBG has been developped and extensively tested using **PacBio HiFi** data.
- MetaMDBG will not work on raw **Nanopore** reads, but error rate is improving quickly, it might work on duplex data in the future. Currently, you have to polish the reads first. For that, you can use [VeChat](https://github.com/HaploKit/vechat) (using Nanopore reads only), or [Ratatosk](https://github.com/DecodeGenetics/Ratatosk) (using Nanopore + Illumina short-reads).
  
## Contig information
Contig information, such as whether it is circular or not, are contained in contig headers in the resulting assembly file.
Examples:

```sh
>ctg1_13x_l
ACGTAGCTTATAGCGAGTATCG...
>ctg2_678x_c
ATTATTGATTAGGGCTATGCAT...
>ctg3_14x_rc
AATTCCGGCGGCGTATTATTAC...
```
Headers are composed of 3 fields separated by underscores.
* Field 1: the name of the contig
* Field 2: estimated coverage for this contig (obtained throught read mapping)
* Field 3: can be "l" (linear), "c" (circular) or "rc" (rescued circular)

Long circular contigs are likely to be complete. Rescued circular are likely to be complete, but it is not guarranted so we recommend using validation methods on them.

## Advanced usage
 
```sh
# Set minimizer length to 16 and use only 0.2% of total k-mers for assembly.
./metaMDBG asm ./outputDir reads.fastq.gz -k 16 -d 0.002

# Stop assembly when reaching a k-mer length of 5000 bps.
./metaMDBG asm ./outputDir reads.fastq.gz -m 5000
```

## Generating an assembly graph

After a successful run of metaMDBG, assembly graph (.gfa) can be generated with the following command.
```sh
./metaMDBG gfa assemblyDir k --contigpath --readpath
```

Assembly dir must be a metaMDBG output dir (the one containing the contig file "contigs.fasta.gz"). The k parameter correspond to the level of resolution of the graph: lower k values will produce graph with high connectivity but shorter unitigs, while higher k graphs will be more fragmented but with longer unitigs. The two optional parameters --contigpath and --readpath allow to generate the path of contigs and reads in the graph respectivelly.

First, display the available k values and their corresponding sequence length in bps (those sequence length in bps are equivalent to the k-mer size that would be used in a traditional de-Brujin graph).
```sh
./metaMDBG gfa ./assemblyDir 0
```

Then, choose a k value and produce the graph (optionnaly add parameters --contigpath and/or --readpath).
```sh
./metaMDBG gfa ./assemblyDir 21
```

MetaMDBG will generate the assembly graph in the GFA format in assemblyDir (e.g. "assemblyGraph_k21_4013bps.gfa").

Note 1) Unitig sequences in the gfa file are not polished, they have the same error rate as in the original reads. Note 2) To generate the unitig sequences, a pass on the original reads that generated the assembly is required, if you have moved the original readsets, you will need to edit the file ./assemblyDir/tmp/input.txt with the new paths.

## Low-memory contig polisher
MetaMDBG contig polisher can be used on any set of contigs. You may be interested by this standalone tool if you have memory issues with existing correction software. Note that the correction method is the same as [Racon](https://github.com/isovic/racon).
```sh
./metaMDBG polish contigs tmpDir reads...

Examples:
./metaMDBG polish assembly.fasta.gz ./tmpDir reads.fastq.gz -t 4                            #Basic usage
./metaMDBG polish assembly.fasta.gz ./tmpDir reads_1.fastq.gz reads_2.fastq.gz -t 4         #Multiple read sets
./metaMDBG polish assembly.fasta.gz ./tmpDir reads_1.fastq.gz reads_2.fastq.gz -t 4 -n 20   #Change maximum read coverage used for correction (here 20x)
```
## Results

Assembly quality and performances on three HiFi PacBio metagenomics samples (using 16 cores).

| Sample | Accession | # bases (Gb) | Wall clock time (h) | Peak memory (GB) | >1Mb near-complete circular contigs | Near-complete MAGs | 
| --- | --- | --- | --- | --- | --- | --- | 
| Human Gut | SRR15275213 | 18.5 | 7 | 6 | 34 | 70 | 
| Anaerobic Digester | ERR10905742 | 64.7  | 13 | 7 | 62 | 130 | 
| Sheep rumen | SRR14289618 | 206.4 | 108 | 22 | 266 | 447 | 

Near-complete: ≥95% completeness and ≤5% contamination (assessed by checkM). Binning was performed with metabat2.

## License

metaMDBG is freely available under the [MIT License](https://opensource.org/license/mit-0/).

## Citation

* Gaetan Benoit, Sebastien Raguideau, Robert James, Adam M. Phillippy, Rayan Chikhi and Christopher Quince [High-quality metagenome assembly from long accurate reads with metaMDBG](https://www.nature.com/articles/s41587-023-01983-6), Nature Biotechnology (2023).
