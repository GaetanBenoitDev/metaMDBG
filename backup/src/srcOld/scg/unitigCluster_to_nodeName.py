

import os, sys, argparse
import matplotlib.pyplot as plt



def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("nodeToContigFile", help="")
    parser.add_argument("contigClusterFile", help="")
    
    args = parser.parse_args()

    nodeToContigFilename = args.nodeToContigFile
    contigClusterFilename = args.contigClusterFile
    clusterFilename = nodeToContigFilename.replace("unitigColor", "unitigCluster")

    nodeToContigFile = open(nodeToContigFilename)
    nodeToContigFile.readline()
    
    unitigToNodenames = {}

    for line in nodeToContigFile:
        line = line.rstrip()
        fields = line.split(",")

        nodeName = int(fields[0])
        unitigIndex = int(fields[1])

        if not unitigIndex in unitigToNodenames:
            unitigToNodenames[unitigIndex] = []

        unitigToNodenames[unitigIndex].append(nodeName)


    outputFile = open(clusterFilename, "w")
    outputFile.write("Name,Color\n")

    contigClusterFile = open(contigClusterFilename)
    contigClusterFile.readline()
    
    for line in contigClusterFile:
        line = line.rstrip()
        fields = line.split(",")

        unitigIndex = int(fields[0])
        cluster = int(fields[1])

        for nodeName in unitigToNodenames[unitigIndex]:
            outputFile.write(str(nodeName) + "," + str(cluster) + "\n")

    outputFile.close()

if __name__ == "__main__":
    main(sys.argv[1:])  