MetaMDBG is a fast and low-memory assembler for long and accurate metagenomics reads (e.g. PacBio HiFi, Nanopore r10.4). It is based on the [minimizer de-Brujin graph](https://github.com/ekimb/rust-mdbg) (MDBG), which have been reimplemetend specifically for metagenomics assembly. MetaMDBG combines an efficient multi-k approach in minimizer-space for dealing with uneven species coverages, and a novel abundance-based filtering method for simplifying strain complexity.

The method nanoMDBG for assembling simplex Nanopore reads (R10.4+) is integrated in metaMDBG.

Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

> [!IMPORTANT]
> 31/07/2024:
> MetaMDBG is in version 1.0: it can now handle nanopore R10.4+ data. It is available on bioconda!
> * Added minimizer-space correction step
> * MetaMDBG parameters changed (options now use explicit names)
> * Contig information format changed in final assemblty fasta file (see Readme - "Contig information") 
> * Removed rescued circular step
> * Removed dependencies to samtools and wfmash
> * Fixed missing "time" dependencies in conda recipe
> * Fixed bug in gfa command


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
Usage:  metaMDBG asm {OPTIONS}

 Basic options:
   --out-dir               Output dir for contigs and temporary files
   --in-hifi               PacBio HiFi read filename(s) (separated by space)
   --in-ont                Nanopore R10.4+ read filename(s) (separated by space)
   --threads               Number of cores [1]

# Nanopore assembly
metaMDBG asm --out-dir ./outputDir/ --in-ont reads.fastq.gz --threads 4
# Hifi assembly
metaMDBG asm --out-dir ./outputDir/ --in-hifi reads.fastq.gz --threads 4
# Multiple sample co-assembly
metaMDBG asm --out-dir ./outputDir/ --in-ont reads_A.fastq.gz reads_B.fastq.gz reads_C.fastq.gz --threads 4
```

MetaMDBG will generate polished contigs in outputDir ("contigs.fasta.gz").
  
## Contig information
Contig information, such as whether it is circular or not, are contained in contig headers in the resulting assembly file.
Examples:

```sh
>ctg112 length=7013 coverage=6 circular=yes
ACGTAGCTTATAGCGAGTATCG...
>ctg37 length=1988 coverage=3 circular=no
ATTATTGATTAGGGCTATGCAT...
>ctg82 length=3824 coverage=13 circular=no
AATTCCGGCGGCGTATTATTAC...
```
Headers are composed of several fields seperated by space.
* **ctgID**:    the name of the contig
* **length**:   the length of the contig in bps
* **coverage**: an estimated read coverage for the contig
* **circular**: whether the contig is circular or no

## Advanced usage
 
```sh
```

## Generating an assembly graph

After a successful run of metaMDBG, assembly graph (.gfa) can be generated with the following command.
```sh
metaMDBG gfa --assembly-dir ./assemblyDir/ --k 21 --contigpath --readpath --threads 4
```

Assembly dir must be a metaMDBG output dir (the one containing the contig file "contigs.fasta.gz"). The --k parameter correspond to the level of resolution of the graph: lower k values will produce graph with high connectivity but shorter unitigs, while higher k graphs will be more fragmented but with longer unitigs. The two optional parameters --contigpath and --readpath allow to generate the path of contigs and reads in the graph respectivelly.

First, display the available k values and their corresponding sequence length in bps (those sequence length in bps are equivalent to the k-mer size that would be used in a traditional de-Brujin graph).
```sh
metaMDBG gfa --assembly-dir ./assemblyDir/ --k 0
```

Then, choose a k value and produce the graph (optionnaly add parameters --contigpath and/or --readpath).
```sh
metaMDBG gfa --assembly-dir ./assemblyDir/ --k 21
```

MetaMDBG will generate the assembly graph in the GFA format in assemblyDir (e.g. "assemblyGraph_k21_4013bps.gfa").

Note 1) Unitig sequences in the gfa file are not polished, they have the same error rate as in the original reads. Note 2) To generate the unitig sequences, a pass on the original reads that generated the assembly is required, if you have moved the original readsets, you will need to edit the file ./assemblyDir/tmp/input.txt with the new paths. Note 3) In nanopore mode, the read-path are not very accurate because of the high error rate, we recommend using actual aligner instead, such as graphAligner.

## Low-memory contig polisher
```sh
```

## Results

Assembly quality and performances on three HiFi PacBio metagenomics samples (using 16 cores).

| Sample | Accession | # bases (Gb) | Wall clock time (h) | Peak memory (GB) | >1Mb near-complete circular contigs | Near-complete MAGs | 
| --- | --- | --- | --- | --- | --- | --- | 
| Human Gut | SRR15275213 | 18.5 | 7 | 6 | 34 | 70 | 
| Anaerobic Digester | ERR10905742 | 64.7  | 13 | 7 | 62 | 130 | 
| Sheep rumen | SRR14289618 | 206.4 | 108 | 22 | 266 | 447 | 

Near-complete: ≥90% completeness and ≤5% contamination (assessed by checkM). Binning was performed with metabat2.

## License

metaMDBG is freely available under the [MIT License](https://opensource.org/license/mit-0/).

## Citation

* Gaetan Benoit, Sebastien Raguideau, Robert James, Adam M. Phillippy, Rayan Chikhi and Christopher Quince [High-quality metagenome assembly from long accurate reads with metaMDBG](https://www.nature.com/articles/s41587-023-01983-6), Nature Biotechnology (2023).
