
rm -rf /pasteur/appa/scratch/gbenoit/run/overlap_test_562/
#rm ~/appa/run/overlap_test_562/tmp/checkpoints/toBasespace.checkpoint

./bin/metaMDBG asm --out-dir /pasteur/appa/scratch/gbenoit/run/overlap_test_562/ --in-ont  ~/appa/data/overlap_test/genome_562_50x/simulatedReads_0.fastq.gz --threads 32

conda run -n quast quast /pasteur/appa/scratch/gbenoit/run/overlap_test_562/contigs.fasta.gz -o /pasteur/appa/scratch/gbenoit/run/overlap_test_562/quast/ -r ~/appa/data/genomes/genomes/562/GCF_000005845.2_ASM584v2_genomic.fna

cat /pasteur/appa/scratch/gbenoit/run/overlap_test_562/quast/report.txt

python3 ~/appa/bin/Autocycler-paper/assess_assembly.py -r ~/appa/data/genomes/genomes/562/GCF_000005845.2_ASM584v2_genomic.fna -a /pasteur/appa/scratch/gbenoit/run/overlap_test_562/contigs.fasta.gz --threads 32

