

import os, sys, argparse, shutil


def main(argv):

	parser = argparse.ArgumentParser()

	#/mnt/gpfs/gaetan/run/experiments/rust-mdbg/AD_origin/binning/bin*.fa
	parser.add_argument("contigs", help="")
	parser.add_argument("reads", help="")
	parser.add_argument("outputDir", help="")

	args = parser.parse_args()

	outputDir = args.outputDir
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)

	alignFilename = outputDir + "/align.bam"
	contigFilename = args.contigs

	command = "minimap2 -t 40 -ax map-hifi " + contigFilename + " " + args.reads + " | samtools sort -o " + alignFilename
	print(command)
	os.system(command)

	command = "samtools index " + alignFilename
	print(command)
	os.system(command)

	command = "java -Xmx600G -jar /mnt/gpfs/gaetan/bin/correction/pilon/pilon-1.24.jar  --genome " + contigFilename + " --pacbio " + alignFilename + " --output contigs_corrected --outdir  " + outputDir + " --fix indels"
	print(command)
	os.system(command)

	#os.remove(alignFilename)


if __name__ == "__main__":
    main(sys.argv[1:])  
