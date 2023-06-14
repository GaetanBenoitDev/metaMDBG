MetaMDBG is a fast and low-memory assembler for long and accurate metagenomics reads. Its core data structure is the [minimizer de-Brujin graph](https://github.com/ekimb/rust-mdbg) (MDBG), which have been reimplemetend specifically for metagenomics assembly. MetaMDBG combines an efficient multi-k approach in minimizer-space for dealing with the wide range of abundances found in metagenomes and a novel abundance-based filtering strategy for simplifying strain complexity.

Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

## Installation

### Conda

Choose an installation directory, then copy-paste all the following commands.
This will create a conda environment, named metaMDBG, with all dependencies installed.
After successful installation, an executable named metaMDBG will appear in ./build/bin.

```
git clone https://github.com/GaetanBenoitDev/metaMDBG.git
cd metaMDBG
conda env create -f conda_env.yml
conda activate metaMDBG
conda env config vars set CPATH=${CONDA_PREFIX}/include:${CPATH}
conda deactivate
conda activate metaMDBG
mkdir build
cd build
cmake ..
make -j 3
```

### Building from source

**Prerequisites**
- gcc 9.4+
- cmake 3.10+
- zlib
- openmp
- minimap2 2.24+
- wfmash
- samtools 1.6+ (using htslib)

```
git clone https://github.com/GaetanBenoitDev/metaMDBG.git
cd metaMDBG
mkdir build
cd build
cmake ..
make -j 3
```

## Usage

```
./metaMDBG asm outputDir reads... {OPTIONS}

	outputDir     Output dir for contigs and temporary files
	reads...      Read filename(s) (separated by space)
	-t            Number of cores [3]
	
Examples:
./metaMDBG asm ./path/to/assemblyDir reads.fastq.gz -t 4                                        #single-sample assembly
./metaMDBG asm ./path/to/assemblyDir reads_A.fastq.gz reads_B.fastq.gz reads_C.fastq.gz -t 4    #co-assembly
```

MetaMDBG will generate polished contigs in outputDir ("contigs.fasta.gz").

## Advanced usage
 
```
# Set minimizer length to 16 and use only 0.2% of total k-mers for assembly.
./metaMDBG asm ./outputDir reads.fastq.gz -k 16 -d 0.002

# Stop assembly when reaching a k-mer length of 5000 bps.
./metaMDBG asm ./outputDir reads.fastq.gz -m 5000
```

## Generating an assembly graph

After a successful run of metaMDBG, assembly graph (.gfa) can be generated with the following command.
```
./metaMDBG gfa assemblyDir k --contigpath --readpath
```

Assembly dir must be a metaMDBG output dir (the one containing the contig file "contigs.fasta.gz"). The k parameter correspond to the level of resolution of the graph: lower k values will produce graph with high connectivity but shorter unitigs, while higher k graphs will be more fragmented but with longer unitigs. The two optional parameters --contigpath and --readpath allow to generate the path of contigs and reads in the graph respectivelly.

First, display the available k values and their corresponding sequence length in bps (those sequence length in bps are equivalent to the k-mer size that would be used in a traditional de-Brujin graph).
```
./metaMDBG gfa ./assemblyDir 0
```

Then, choose a k value and produce the graph (optionnaly add parameters --contigpath and/or --readpath).
```
./metaMDBG gfa ./assemblyDir 21
```

MetaMDBG will generate the assembly graph in the GFA format in assemblyDir (e.g. "assemblyGraph_k21_4013bps.gfa").

Note 1) Unitig sequences in the gfa file are not polished, they have the same error rate as in the original reads. Note 2) To generate the unitig sequences, a pass on the original reads that generated the assembly is required, if you have moved the original readsets, you will need to edit the file ./assemblyDir/tmp/input.txt with the new paths.

## Low-memory contig polisher
MetaMDBG contig polisher can be used on any set of contigs. You may be interested by this standalone tool if you have memory issues with exsting correction software. Note that the correction method is the same as [Racon](https://github.com/isovic/racon).
```
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

Near-complete: ≥95% completeness and ≤5 contamination (assessed by checkM)
MAGs: computed by metabat2

## License

metaMDBG is freely available under the [MIT License](https://opensource.org/license/mit-0/).

## Citation

