# pyabpoa: abPOA Python interface
## Introduction
pyabpoa provides an easy-to-use interface to [abPOA](https://github.com/yangao07/abPOA), it contains all the APIs that can be used to perform MSA for a set of sequences and consensus calling from the final alignment graph.

## Installation

### Install pyabpoa with pip

pyabpoa can be installed with pip:

```
pip install pyabpoa
```

### Install pyabpoa from source
Alternatively, you can install pyabpoa from source (cython is required):
```
git clone --recursive https://github.com/yangao07/abPOA.git
cd abPOA
make install_py
```

## Examples
The following code illustrates how to use pyabpoa.
```
import pyabpoa as pa
a = pa.msa_aligner()
seqs=[
'CCGAAGA',
'CCGAACTCGA',
'CCCGGAAGA',
'CCGAAGA'
]
res=a.msa(seqs, out_cons=True, out_msa=True, out_pog='pog.png', incr_fn='') # perform multiple sequence alignment 
                                                                # generate a figure of alignment graph to pog.png

for seq in res.cons_seq:
    print(seq)  # print consensus sequence

res.print_msa() # print row-column multiple sequence alignment in PIR format
```
You can also try the example script provided in the source folder:
```
python ./python/example.py
```


## APIs

### Class pyabpoa.msa_aligner
```
pyabpoa.msa_aligner(aln_mode='g', ...)
```
This constructs a multiple sequence alignment handler of pyabpoa, it accepts the following arguments:

* **aln_mode**: alignment mode. 'g': global, 'l': local, 'e': extension; default: **'g'**
* **is_aa**: input is amino acid sequence; default: **False**
* **match**: match score; default: **2**
* **mismatch**: match penaty; default: **4**
* **score_matrix**: scoring matrix file, **match** and **mismatch** are not used when **score_matrix** is used; default: **''**
* **gap_open1**: first gap opening penalty; default: **4**
* **gap_ext1**: first gap extension penalty; default: **2**
* **gap_open2**: second gap opening penalty; default: **24**
* **gap_ext2**: second gap extension penalty; default: **1**
* **extra_b**: first adaptive banding paremeter; set as < 0 to disable adaptive banded DP; default: **10**
* **extra_f**: second adaptive banding paremete; the number of extra bases added on both sites of the band is *b+f\*L*, where *L* is the length of the aligned sequence; default : **0.01**

The `msa_aligner` handler provides one method which performs multiple sequence alignment and takes four arguments:
```
pyabpoa.msa_aligner.msa(seqs, out_cons, out_msa, out_pog='', incr_fn='')
```

* **seqs**: a list variable containing a set of input sequences; **positional**
* **out_cons**: a bool variable to ask pyabpoa to generate consensus sequence; **positional**
* **out_msa**: a bool variable to ask pyabpoa to generate RC-MSA; **positional**
* **max_n_cons**: maximum number of consensus sequence to generate; default: **1**
* **min_freq**: minimum frequency of each consensus to output (effective when **max_n_cons** > 1); default: **0.3**
* **out_pog**: name of a file (`.png` or `.pdf`) to store the plot of the final alignment graph; default: **''**
* **incr_fn**: name of an existing graph (GFA) or MSA (FASTA) file, incrementally align sequence to this graph/MSA; default: **''**

### Class pyabpoa.msa_result
```
pyabpoa.msa_result(seq_n, cons_n, cons_len, ...)
```
This class describes the information of the generated consensus sequence and the RC-MSA. The returned result of `pyabpoa.msa_aligner.msa()` is an object of this class that has the following properties:

* **n_seq**: number of input aligned sequences
* **n_cons**: number of generated consensus sequences (generally 1, could be 2 or more if **max_n_cons** is set as > 1)
* **clu_n_seq**: an array of sequence cluster size
* **cons_len**: an array of consensus sequence length(s)
* **cons_seq**: an array of consensus sequence(s)
* **cons_cov**: an array of consensus sequence coverage for each base
* **msa_len**: size of each row in the RC-MSA
* **msa_seq**: an array containing `n_seq`+`n_cons` strings that demonstrates the RC-MSA, each consisting of one input sequence and several `-` indicating the alignment gaps. 

`pyabpoa.msa_result()` has a function of `print_msa` which prints the RC-MSA to screen.

```
pyabpoa.msa_result().print_msa()
```
