
rm -rf /pasteur/appa/scratch/gbenoit/run/overlap_test_562/

./bin/metaMDBG asm --out-dir /pasteur/appa/scratch/gbenoit/run/overlap_test_562/ --in-ont ~/appa/data/overlap_test/genome_562_50x/simulatedReads_0.fastq.gz --threads 32 

conda run -n quast quast /pasteur/appa/scratch/gbenoit/run/overlap_test_562/contigs.fasta.gz -o /pasteur/appa/scratch/gbenoit/run/overlap_test_562/quast/ -r ~/appa/data/genomes/genomes/562/GCF_000005845.2_ASM584v2_genomic.fna

cat /pasteur/appa/scratch/gbenoit/run/overlap_test_562/quast/report.txt
