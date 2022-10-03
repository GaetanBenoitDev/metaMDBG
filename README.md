metaMDBG is a lightweight assembler for long and accurate metagenomics reads.  
Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

## Dependencies

```
- gcc 9.4+
- cmake 3.10+
- zlib
- openmp
- minimap2
```

## Installation

```
git clone https://github.com/GaetanBenoitDev/metaMDBG.git
cd metaMDBG
mkdir build
cd build
cmake ..
make -j 3
```

After successful installation, an executable named metaMDBG will appear in build/bin.


## Usage

```
./metaMDBG asm outputDir reads... {OPTIONS}

	outputDir     Output dir for contigs and temporary files
	reads...      Read filename(s) (separated by space)
	-t            Number of cores [3]
```

MetaMDBG will generate polished contigs in outputDir ("contigs_polished.fasta.gz").

## Purging strain duplication (experimental)

