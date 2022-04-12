




import os, sys, argparse


mdbgFilename = "/mnt/gpfs/gaetan/tools/MdbgAssembler/build/bin/mdbgAsmMeta"
drepScriptFilename = "/mnt/gpfs/gaetan/scripts/annotation/annotation2/drepStats.py"


def main(argv):
    
    parser = argparse.ArgumentParser()

    #/mnt/gpfs/gaetan/run/experiments/rust-mdbg/AD_origin/binning/bin*.fa
    parser.add_argument("ilong", help="input HiFi")
    parser.add_argument("ishort", help="input short reads")
    parser.add_argument("asmDir", help="dir for assembly files")
    parser.add_argument("binDir", help="output dir for bins")
    parser.add_argument("t", help="nb threads")
    #parser.add_argument("drep", help="path to drep script")
    
    args = parser.parse_args()
    inputFilenameLong = args.ilong
    inputFilenameShort = args.ishort
    nbCores = args.t
    asmDir = args.asmDir
    binDir = args.outDir

    contigsFilename = asmDir + "/contigs_81.fasta.gz"
    binRegex = "\"" + binDir + "/bin_*.fasta" + "\""
    drepDir = binDir + "/drep"
    kmerCoverageFilename = asmDir + "/kmerCoverages_k81.tsv"

    command = mdbgFilename + " asm -i " + inputFilenameLong + " -o " + asmDir + " -t " + str(nbCores)
    command = mdbgFilename + " countKmer -o " + kmerCoverageFilename + " -i " + inputFilenameShort + " -c " + contigsFilename + " -k 61 -t " + str(nbCores)
    command = mdbgFilename + " bin " + contigsFilename + " " + asmDir + " " + binDir  + " -a " + kmerCoverageFilename
    command = "python3 " + drepScriptFilename + " " + binRegex + " " + drepDir


if __name__ == "__main__":
    main(sys.argv[1:])  
