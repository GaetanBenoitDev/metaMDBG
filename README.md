MetaMDBG is a fast, memory-efficient assembler designed for long and accurate metagenomics reads (e.g. PacBio HiFi, Nanopore R10). It is optimized for metagenomes, but also works well on bacterial isolates samples. Up-to-date benchmarks are provided in the Results section below.

The method [nanoMDBG](https://www.nature.com/articles/s41467-026-69760-y) for assembling Nanopore R10 simplex reads is integrated in metaMDBG.

Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

## News
Feb 2026:
MetaMDBG v1.3 is out!
* Fixed clipping events
* Fixed zero-coverage regions
* Fixed chimeric contigs
* Improved performances
  
Check out the new results below.


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
conda activate metamdbg1.4
conda env config vars set CPATH=${CONDA_PREFIX}/include:${CPATH}
conda deactivate

# Activate metaMDBG environment
conda activate metamdbg1.4

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
  
After successful installation, an executable named metaMDBG will appear in ./build/bin.

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
Contig information are contained in contig headers in the resulting fasta assembly file.
Example:

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

## Resume an existing run (checkpoint system)

If an assembly run stops for any reason, simply resubmit the same command.
MetaMDBG will automatically skip completed steps and resume from the last checkpoint.

## Advanced usage
 
```sh
# Filter out reads with low average per-base quality (using phred score)
metaMDBG asm --out-dir ./outputDir/ --in-ont reads.fastq.gz --min-read-quality 10

# Filter output contigs (by length and by coverage)
metaMDBG asm --out-dir ./outputDir/ --in-ont reads.fastq.gz --min-contig-length 500 --min-contig-coverage 2

# Skip correction step (useful if using corrected reads)
metaMDBG asm --out-dir ./outputDir/ --in-ont reads.fastq.gz --skip-correction

# Filter out unique k-min-mers to improve performances.
# Useful for scaling to very large datasets, but may reduce assembly quality and completeness.
# By default, metaMDBG attempts to rescue low-abundance genomic k-min-mers.
metaMDBG asm --out-dir ./outputDir/ --in-ont reads.fastq.gz --min-abundance 2

# Stop assembly after reaching k-th iteration.
metaMDBG asm --out-dir ./outputDir/ --in-ont reads.fastq.gz --max-k 11
```

## Generating an assembly graph

<details><summary>See details</summary>

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

</details>

## Results

[<img src="https://github.com/GaetanBenoitDev/metaMDBG/blob/main/results/fig_metaMDBG_1.3.jpg" width="75%" />](https://github.com/GaetanBenoitDev/metaMDBG/blob/main/results/fig_metaMDBG_1.3.jpg)

Source data: [mags.tsv](https://github.com/GaetanBenoitDev/metaMDBG/blob/main/results/mags.tsv) [errors.tsv](https://github.com/GaetanBenoitDev/metaMDBG/blob/main/results/errors.tsv) [perf.tsv](https://github.com/GaetanBenoitDev/metaMDBG/blob/main/results/perf.tsv)

Alignment and binning were performed with minimap2 and SemiBin2. Completeness and contamination were measured with checkM2 (near-complete: ≥90% completeness and ≤5% contamination, Medium: ≥50% completeness and ≤5% contamination). Clipping events and zero-coverage regions were identified using the anvi-script-find-misassembly program from the Anvi’o platform. All assemblers were run with 32 cores. 


| Sample | Accession | # bases (Gb) | N50 read length (kb) | Average quality score |  
| --- | --- | --- | --- | --- | 
| Human Gut 1 (ONT) | ERR15285694 | 50 | 7.8 | 23.2 | 
| Human Gut 2 (ONT) | SRR29980972 | 77 | 27.2 | 17.3 | 
| Oral (ONT) | DRR582205 | 24 | 15 | 21.7 | 
| Soil Microflora (ONT) | ERR11523665 | 103 | 5.4 | 17.1 | 
| Human Gut 1 (HiFi) | ERR15289675 | 50 | 8.9 | 34 | 
| Human Gut 2 (HiFi) | SRR15275213 | 18.5 | 11.4 | 45 | 
| Anaerobic Digester (HiFi) | ERR10905743 | 67 | 10.2 | 40.6 | 
| Sea Water (HiFi) | ERR9769281 | 22 | 8.2 | 35 | 

## Results on bacterial isolates

| Sample | Average quality score | # reference genomes | # contigs | Genome fraction (%) | Duplication ratio | # mismatches per 100 kbp | # indels per 100 kbp | # misassemblies | # local misassemblies | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| [Escherichia_coli_1](https://journals.asm.org/doi/full/10.1128/mra.00006-26) | 20.97 | 6 | 6 | 99.999 | 1.000 | 0.00 | 0.14 | 0 | 0 | 
| [Mycobacterium_avium_1](https://journals.asm.org/doi/full/10.1128/mra.00014-26) | 28.08 | 1 | 1 | 100.000 | 1.000 | 0.00 | 0.05 | 0 | 0 | 
| [Phascolarctobacterium_faecium_1](https://journals.asm.org/doi/full/10.1128/mra.00789-25) | 15.45 | 1 | 1 | 100.000 | 1.000 | 8.60 | 12.90 | 0 | 0 | 
| [Escherichia_coli_2](https://journals.asm.org/doi/full/10.1128/mra.00746-25) | 16.23 | 5 | 4 | 99.905 | 1.000 | 1.09 | 1.47 | 0 | 1 | 
| [Saccharolobus_islandicus_1](https://journals.asm.org/doi/full/10.1128/mra.00626-25) | 27.54 | 1 | 1 | 99.838 | 1.000 | 0.00 | 0.00 | 0 | 0 | 
| [Agarivorans_1](https://journals.asm.org/doi/full/10.1128/mra.00577-25) | 20.73 | 1 | 1 | 100.000 | 1.000 | 0.28 | 0.12 | 0 | 0 | 


MetaMDBG was applied to bacterial isolate samples (each containing a single bacterial genome and, in some cases, multiple plasmids). The linked studies provide the original sequencing reads and the reference genomes assembled using Autocycler. Assembly quality was evaluated using QUAST. 

## License

metaMDBG is freely available under the [MIT License](https://opensource.org/license/mit-0/).

## Citation

* Gaetan Benoit, Sebastien Raguideau, Robert James, Adam M. Phillippy, Rayan Chikhi and Christopher Quince [High-quality metagenome assembly from long accurate reads with metaMDBG](https://www.nature.com/articles/s41587-023-01983-6), Nature Biotechnology (2023).

* Gaetan Benoit, Robert James, Sebastien Raguideau, Georgina Alabone, Tim Goodall, Rayan Chikhi and Christopher Quince [High-quality metagenome assembly from nanopore reads with nanoMDBG](https://www.nature.com/articles/s41467-026-69760-y), Nature Communications (2026).
