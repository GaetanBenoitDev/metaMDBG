metaMDBG is a lightweight assembler for metagenomics long accurate reads.

Developper: Gaëtan Benoit

Contact: gaetanbenoitdev at gmail dot com

## Dependencies

```
- gcc 9.4+
- cmake 3.10+
- zlib
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
	reads...      Input filename(s) (separated by space)
	-l            Minimizer length [13]
	-d            Minimizer density [0.005]
	-t            Number of cores [3]
	--nofilter    Disable unique kminmer filter prior to graph construction
```
