
import os, sys, argparse, glob
import matplotlib.pyplot as plt
import numpy as np


def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", help="input stat dir")
    
    args = parser.parse_args()

    execute(args.dir)

def execute(dirpath):

    for filename in glob.glob(os.path.join(dirpath, "*.txt")):
        outputFilename = os.path.splitext(filename)[0]
        outputFilename += ".png"
        print(outputFilename)

        abundances = []

        for line in open(filename):
            line = line.rstrip()
            fields = line.split(" ")
            for field in fields:
                abundances.append(int(field))

        median = np.median(abundances)
        print(median)

        abundances_cleaned = []
        for abundance in abundances:
            if(abundance > median*2): continue
            abundances_cleaned.append(abundance)

        nbBins = int(median*2) #int(median+median+0.75 - median-median*0.75)
        #print(nbBins)

        plt.clf()
        #print(abundances)
        plt.hist(abundances_cleaned, density=True, bins=nbBins)  # density=False would make counts
        plt.ylabel('Probability')
        plt.xlabel('Data')
        #plt.xlim([median-median*0.75, median+median+0.75])
        plt.savefig(outputFilename)



if __name__ == "__main__":
    main(sys.argv[1:])  
