




import os, sys, argparse





def main(argv):
    
    parser = argparse.ArgumentParser()

    #/mnt/gpfs/gaetan/run/experiments/rust-mdbg/AD_origin/binning/bin*.fa
    parser.add_argument("i", help="input filename")
    parser.add_argument("tmpDir", help="dir for temporary files")
    parser.add_argument("t", help="nb threads")
    
    args = parser.parse_args()
    inputFilename = args.i
    nbCores = args.t
    tmpDir = args.tmpDir

    command = "nohup ./bin/mdbgAsmMeta asm -i " + inputFilename + " -o " + tmpDir + " -t " + str(nbCores)
    command = "mdbgAsmMeta countKmer -o ../../../run/mdbg/AD_origin_noSelfCycle/kmerCoverages_k81.tsv -i ../../../data/HiFi_AD/input_shortreads_all.txt -c ../../../run/mdbg/AD_origin_noSelfCycle/contigs_81.fasta.gz -k 61 -t 64"



if __name__ == "__main__":
    main(sys.argv[1:])  
