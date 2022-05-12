import pyabpoa as pa

#@parameters of msa_aligner:
#   aln_mode='g' # g: global, l: local, e: extension
#   is_aa=False # set as True if input is amino acid sequence
#   score_matrix='' # file of score matrix, e.g. HOXD70.mtx/BLOSUM62.mtx
#   match=2
#   mismatch=4
#   gap_open1=4
#   gap_open2=24
#   gap_ext1=2
#   gap_ext2=1
#   extra_b = 10 # 1st part of extra band, -1 to disable banded DP
#   extra_f = 0.01  # 2nd part of eatra band. w = extra_b+extra_f*L (L is sequence length)
#   max_n_cons=1 # to output at most N cons, set max_n_cons as N
#   min_freq=0.3 # minimum frequence of each consensus to output for diploid data

# construct msa aligner
a = pa.msa_aligner()

print("==== First exmaple: 2 consensus sequences ====\n")
# for multiple consensus
seqs=[
 'CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT',
 'CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT'
 ]

#@parameters of msa
#seqs: multiple sequences
out_cons=True # generate consensus sequence, set as False to disable
out_msa=True # generate row-column multiple sequence alignment, set as False to disable
out_pog="example1.png" # generate plot of alignment graph, set None to disable
max_n_cons = 2

# multiple sequence alignment for 'seqs'
res=a.msa(seqs, out_cons=out_cons, out_msa=out_msa, max_n_cons=max_n_cons, out_pog=out_pog)

# output result
if out_cons:
    for i in range(res.n_cons):
        print(">Consensus_sequence_{}".format(i+1))
        print(res.cons_seq[i])
if out_msa:
    res.print_msa()


print("\n\n==== Second exmaple: 1 consensus sequence ====\n")
seqs=[
'CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA',
'CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC',
'CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC',
'CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC',
'CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC',
'CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC',
'CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC',
'CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC',
'CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC',
'CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT'
]

#@parameters of msa
#seqs: multiple sequences
out_cons=True # generate consensus sequence, set as False to disable
out_msa=True # generate row-column multiple sequence alignment, set as False to disable
out_pog="example2.png" # generate plot of alignment graph, set None to disable
max_n_cons = 1

# multiple sequence alignment for 'seqs'
res=a.msa(seqs, out_cons=out_cons, out_msa=out_msa, max_n_cons=max_n_cons, out_pog=out_pog)

# output result
if out_cons:
    for i in range(res.n_cons):
        print(">Consensus_sequence_{}".format(i+1))
        print(res.cons_seq[i])
if out_msa:
    for i in range(res.n_seq):
        print(">Seq_{}".format(i+1))
        print(res.msa_seq[i])
    for i in range(res.n_cons):
        print(">Consensus_sequence_{}".format(i+1))
        print(res.msa_seq[res.n_seq+i])