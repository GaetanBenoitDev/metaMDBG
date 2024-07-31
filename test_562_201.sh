
#/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/run/overlap_test_201/
./bin/metaMDBG asm /pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/run/overlap_test_201/ /pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz /pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/data/overlap_test/genome_562_50x/simulatedReads_0.fastq.gz -t 16

conda run -n quast quast /pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/run/overlap_test_201/contigs.fasta.gz -o /pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/run/overlap_test_201/quast/ -r /pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/data/genomes/genomes/201/GCF_000816365.1_ASM81636v1_genomic.fna

