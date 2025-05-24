
./bin/metaMDBG asm --out-dir /pasteur/appa/scratch/gbenoit/run/overlap_test_201/ --in-hifi ~/appa/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz --threads 8 --kmer-size 15

conda run -n quast quast /pasteur/appa/scratch/gbenoit/run/overlap_test_201/contigs.fasta.gz -o /pasteur/appa/scratch/gbenoit/run/overlap_test_201/quast/ -r /pasteur/appa/scratch/gbenoit/data/genomes/genomes/201/GCF_000816365.1_ASM81636v1_genomic.fna

cat /pasteur/appa/scratch/gbenoit/run/overlap_test_201/quast/report.txt
