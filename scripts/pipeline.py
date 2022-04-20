




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
    parser.add_argument("kminmerLength", help="kminmer length")
    parser.add_argument("t", help="nb threads")
    #parser.add_argument("drep", help="path to drep script")
    
    args = parser.parse_args()
    inputFilenameLong = args.ilong
    inputFilenameShort = args.ishort
    nbCores = args.t
    asmDir = args.asmDir
    binDir = args.binDir

    contigsFilename = asmDir + "/contigs_" + args.kminmerLength + ".fasta.gz"
    binRegex = "\"" + binDir + "/bins_/bin_*.fasta" + "\""
    drepDir = binDir + "/drep"
    kmerCoverageFilename = asmDir + "/kmerCoverages_k" + args.kminmerLength + ".tsv"

    command = mdbgFilename + " asm -i " + inputFilenameLong + " -o " + asmDir + " -t " + str(nbCores)
    #execute_command(command)

    command = mdbgFilename + " countKmer -o " + kmerCoverageFilename + " -i " + inputFilenameShort + " -c " + contigsFilename + " -k 61 -t " + str(nbCores)
    execute_command(command)

    command = mdbgFilename + " bin " + contigsFilename + " " + asmDir + " " + binDir  + " -a " + kmerCoverageFilename
    execute_command(command)
    
    command = "python3 " + drepScriptFilename + " " + binRegex + " " + drepDir
    execute_command(command)

def execute_command(command):

    print(command)
    ret = os.system(command)

    if ret != 0:
        print("Error: ", command)
        exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])  
