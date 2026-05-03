
rm -rf /pasteur/appa/scratch/gbenoit/run/overlap_test_201/
#rm ~/appa/run/overlap_test_201/tmp/checkpoints/toBasespace.checkpoint
#cp -r /pasteur/appa/scratch/gbenoit/run/overlap_test_201_1.3/ /pasteur/appa/scratch/gbenoit/run/overlap_test_201/

#rm ~/appa/run/overlap_test_201/tmp/checkpoints/k*
#rm ~/appa/run/overlap_test_201/tmp/checkpoints/derep*
#rm ~/appa/run/overlap_test_201/tmp/checkpoints/toBasespace*
#rm -rf ~/appa/run/overlap_test_201/tmp/pass_k*
#rm -rf ~/appa/run/overlap_test_201/tmp/smallContigs/
#rm ~/appa/run/overlap_test_201/tmp/unitig*
#rm ~/appa/run/overlap_test_201/tmp/contig*

./bin/metaMDBG asm --out-dir /pasteur/appa/scratch/gbenoit/run/overlap_test_201/ --in-ont ~/appa/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz --threads 32 --all-assembly-graph

conda run -n quast quast /pasteur/appa/scratch/gbenoit/run/overlap_test_201/contigs.fasta.gz -o /pasteur/appa/scratch/gbenoit/run/overlap_test_201/quast/ -r /pasteur/appa/scratch/gbenoit/data/genomes/genomes/201/GCF_000816365.1_ASM81636v1_genomic.fna

cat /pasteur/appa/scratch/gbenoit/run/overlap_test_201/quast/report.txt
