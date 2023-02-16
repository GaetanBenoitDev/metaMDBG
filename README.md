MetaMDBG is a lightweight assembler for long and accurate metagenomics reads.

Developper: Gaëtan Benoit  
Contact: gaetanbenoitdev at gmail dot com

## Dependencies

```
- gcc 9.4+
- cmake 3.10+
- zlib
- openmp
- minimap2 2.24+
- wfmash
- samtools 1.6+
```

## Conda Installation

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

After successful installation, an executable named metaMDBG will appear in build/bin.

## Installation from source
Required dependencies to be installed manually.

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
```

MetaMDBG will generate contigs in outputDir ("contigs.fasta.gz").


## License

rust-mdbg is freely available under the [MIT License](https://opensource.org/license/mit-0/).

## Citation

ongoing