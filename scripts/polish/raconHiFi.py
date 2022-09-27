

import os, sys, argparse, shutil


def main(argv):

	parser = argparse.ArgumentParser()

	#/mnt/gpfs/gaetan/run/experiments/rust-mdbg/AD_origin/binning/bin*.fa
	parser.add_argument("contigs", help="Contig filename to be corrected")
	parser.add_argument("reads", help="Read filename")
	parser.add_argument("outputFilename", help="Output contig filename")
	parser.add_argument("nbCores", help="")

	args = parser.parse_args()

	nbCores = int(args.nbCores)
	alignFilename = args.contigs + "__.sam.gz"

	outputFilename = args.outputFilename
	if outputFilename == args.contigs: sys.exit(1)

	outputGz = ".gz" in outputFilename

	command = "minimap2 -ax map-hifi " + args.contigs + " " + args.reads + " -t " + str(nbCores) + " | gzip -c - > " + alignFilename
	print(command)
	os.system(command)

	outputCommand = " > " + outputFilename
	if outputGz: outputCommand = " | gzip -c - > " + outputFilename


	command = "racon -t " + str(nbCores) + " " + args.reads + " " + alignFilename + " " + args.contigs + " " + outputCommand #+ os.path.splitext(args.contigs)[0] + "_corrected.fasta"
	print(command)
	os.system(command)

	os.remove(alignFilename)


if __name__ == "__main__":
    main(sys.argv[1:])  

