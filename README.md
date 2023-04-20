MetaMDBG is a lightweight assembler for long and accurate metagenomics reads.

Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

## Installation

### Conda

Choose a directory to install metaMDBG, then copy-paste all the following commands.
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

### Basic usage
```
./metaMDBG asm outputDir reads... {OPTIONS}

	outputDir     Output dir for contigs and temporary files
	reads...      Read filename(s) (separated by space)
	-t            Number of cores [3]
```

MetaMDBG will generate polished contigs in outputDir ("contigs.fasta.gz").

### Advanced usage
 
- Set minimizer length to 16 and use 0.2% of k-mer for assembly.
```
./metaMDBG asm ./outputDir reads.fastq.gz -k 16 -d 0.002
```
- Set minimizer length to 16 and use 0.2% of k-mer for assembly.
```
./metaMDBG asm ./outputDir reads.fastq.gz -k 16 -d 0.002
```

## License

metaMDBG is freely available under the [MIT License](https://opensource.org/license/mit-0/).

## Citation

ongoing
