




import os, sys, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import marker_gene_utils
import statistics





def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("contig", help="")
    parser.add_argument("nodeToContigFile", help="")
    parser.add_argument('--isSingleSpecies', dest='isSingleSpecies', action='store_true')


    #parser.add_argument("out", help="Output filename")
    
    #parser.add_argument("csv", help="output unitig coverage file (.csv)")
    
    args = parser.parse_args()

    contigFilename = args.contig
    nodeToContigFilename = args.nodeToContigFile
    isSingleSpecies = args.isSingleSpecies

    nbCores = 4
    mg_threshold = 0.5 #length threshold to consider marker genes
    min_length = 1000 #minimum length of contigs to consider for compositional probability

    contigLengths = {}
    #contigIndex = 0
    contigName_to_contigIndex = {}
    for header, seq in SimpleFastaParser(open(contigFilename)):
        contigIndex = int(header.replace("ctg", ""))
        contigName_to_contigIndex[header] = contigIndex
        #contigIndex += 1
        contigLengths[contigIndex] = len(seq)

    print(contigName_to_contigIndex)
    print(contigLengths.keys())

    scgClusterFile_contig = open(contigFilename + ".scgCluster.txt", "w")
    scgClusterFile_contig.close()

    scgSequencesFilename = contigFilename + ".scgSequences.ffn"
    if os.path.exists(scgSequencesFilename):
        os.remove(scgSequencesFilename)

    #print(isSingleSpecies)

    if ".gz" in contigFilename:
        print("Besoin de decomrpesser le contig.fasta.gz")
        exit(1)

    if isSingleSpecies:
        checkIsSingleSpecies(contigFilename, nodeToContigFilename, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths)
    else:
        annotate(contigFilename, nodeToContigFilename, False, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths)
        computeScgClusters(contigFilename, nodeToContigFilename, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths)

def getUnitigToNodenames(nodeToContigFilename):

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

    return unitigToNodenames

def checkIsSingleSpecies(contigFilename, nodeToContigFilename, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths):

    resultFile = open(contigFilename + "_isSingleSpecies.txt", "w")

    splits_fasta = []
    nbUnitigs = 500

    nbUnitigSeen = 0

    contigFilenameTmp = contigFilename + ".part"
    fastaFile = open(contigFilenameTmp, "w")
    marker_contig_counts_all = {}

    scgSequencesFilename = contigFilenameTmp + ".scgSequences.ffn"
    if os.path.exists(scgSequencesFilename):
        os.remove(scgSequencesFilename)

    for header, seq in SimpleFastaParser(open(contigFilename)):

        #print(nbUnitigSeen, len(seq))

        fastaFile.write(">" + header + "\n")
        fastaFile.write(seq + "\n")
        nbUnitigSeen += 1

        if nbUnitigSeen > nbUnitigs:
            fastaFile.close()
            annotatePass(contigFilenameTmp, nodeToContigFilename, True, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths, resultFile, False)

            os.remove(contigFilenameTmp)
            fastaFile = open(contigFilenameTmp, "w")
            nbUnitigSeen = 0

    fastaFile.close()
    if nbUnitigSeen > 0:
        annotatePass(contigFilenameTmp, nodeToContigFilename, True, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths, resultFile, True)

    os.remove(contigFilenameTmp)

    print(marker_contig_counts_all)

    resultFile.write("1")
    resultFile.close()


def annotatePass(contigFilenameTmp, nodeToContigFilename, isAnnotation, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths, resultFile, isLastPass):

    scgCounts = annotate(contigFilenameTmp, nodeToContigFilename, isAnnotation, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths)
        
    completeness, contamination = computeCompletenessContamination2(scgCounts)
    print("Completeness:", completeness)
    print("Contamination:", contamination)
    if not isLastPass and contamination <= 0.05: return
    
    scgClusterCounts = computeScgClusters(contigFilenameTmp, nodeToContigFilename, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths)

    completeness, contamination = computeCompletenessContamination(scgClusterCounts)
    print("Completeness:", completeness)
    print("Contamination:", contamination)

    if contamination > 0.05:
        resultFile.write("0")
        resultFile.close()
        exit(0)

    """
    marker_contig_counts = annotate(contigFilenameTmp, nodeToContigFilename, isAnnotation, mg_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths)
    for scgName, count in marker_contig_counts.items():
        if not scgName in marker_contig_counts_all:
            marker_contig_counts_all[scgName] = 0
        marker_contig_counts_all[scgName] += count

    #print(marker_contig_counts_all)
    print(len(marker_contig_counts_all))

    nbSpecies = 0
    if len(marker_contig_counts_all) > 30:
        nbSpecies = int(statistics.median(marker_contig_counts_all.values()))
    print("Predicted nb species: ", nbSpecies)

    if nbSpecies != 1:
        resultFile.write("0")
        resultFile.close()
        exit(0)
    """


def annotate(contigFilename, nodeToContigFilename, annotateOnly, mg_length_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths):
    





    scriptDir = os.path.dirname(os.path.realpath(__file__))

    # Run FragGeneScan and HMMER if .hmmout file is not present
    marker_gene_utils.scan_for_marker_genes(contigs_file=contigFilename, nthreads=nbCores, markerURL=scriptDir+"/marker.hmm")


    unitigToNodenames = getUnitigToNodenames(nodeToContigFilename)

    # Get contigs with single-copy marker genes and count of contigs for each single-copy marker gene
    marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes(
        contigs_file=contigFilename,
        contig_names_rev=contigName_to_contigIndex,
        mg_length_threshold=mg_length_threshold,
        contig_lengths=contigLengths,
        min_length=min_length)

    #all_contig_markers = marker_gene_utils.get_all_contigs_with_marker_genes(
    #    contigs_file=contigFilename,
    #    contig_names_rev=contigName_to_contigIndex,
    #    mg_length_threshold=mg_threshold)

    print(marker_contigs)
    print("----")
    print(marker_contig_counts)
    print("----")
    print(contig_markers)

    #if annotateOnly: return marker_contig_counts

    scgIndex = 0
    scgName_to_index = {}
    for scgName in marker_contigs.keys():
        scgName_to_index[scgName] = scgIndex
        scgIndex += 1

    #print(scgName_to_index)

    scfFile = open(nodeToContigFilename.replace("Color", "SCG"), "w")
    scfFile.write("Name,Color\n")

    nodenamesUsed = {}
    for unitigName, scgNames in contig_markers.items():
        for scgName in scgNames:
            for nodeName in unitigToNodenames[unitigName]:
                if nodeName in nodenamesUsed: continue

                scfFile.write(str(nodeName) + "," + str(scgName_to_index[scgName]) + "\n")
                nodenamesUsed[nodeName] = True
                break

    scfFile.close()



    scgSequencesNames = {}

    hmmFilename = contigFilename + ".hmmout"
    for line in open(hmmFilename):
        line = line.rstrip()

        if line[0] == "#": continue

        fields = line.split()
        #print(fields)
        scgSequenceName = fields[0]
        cogName = fields[3]
        cogName = cogName.replace("_", "")


        #parsercontig = strings[0]

        # Marker gene name
        #marker_gene = strings[3]

        # Marker gene length
        marker_gene_length = int(fields[5])

        # Mapped marker gene length
        mapped_marker_length = int(fields[16]) - int(fields[15])

        name_strings = scgSequenceName.split("_")
        name_strings = name_strings[:len(name_strings)-3]

        # Contig name
        contig_name = "_".join(name_strings)

        contig_num = contigName_to_contigIndex[contig_name]
        contig_length = contigLengths[contig_num]

        if contig_length >= min_length and mapped_marker_length > marker_gene_length*mg_length_threshold:

            scgSequencesNames[scgSequenceName] = cogName

    scgSequencesFilename = contigFilename + ".scgSequences.ffn"
    scgSequenceFile = open(scgSequencesFilename, "a")

    geneFilename = contigFilename + ".frag.ffn"
    for header, seq in SimpleFastaParser(open(geneFilename)):
        if header in scgSequencesNames:
            
            cogName = scgSequencesNames[header]
            header += "_" + cogName
            scgSequenceFile.write(">" + header + "\n")
            scgSequenceFile.write(seq + "\n")

    #print(len(scgSequencesNames))
    scgSequenceFile.close()


    """
    nbSpecies = 0
    if len(marker_contig_counts) > 0:
        nbSpecies = int(statistics.median(marker_contig_counts.values()))

    print("Predicted nb species: ", nbSpecies)
    if nbSpecies > 1:
        
        distanceMatrixFilename = contigFilename + "._distanceMatrixComposition.csv.gz";
        compositionCommand = "/home/gats/workspace/tools/computeUnitigAbundance/build/bin/Bloocoo -file " + contigFilename + " -abd ~/workspace/run/shortreads/unitigs/depth_input.txt"
        print(compositionCommand)
        ret = os.system(compositionCommand)
        if ret != 0:
            print("Command failed: ", ret)

        contigClusterFilename = distanceMatrixFilename + ".csv"
        clusterCommand = "Rscript " + scriptDir + "/clusterContigsComposition.r " + distanceMatrixFilename + " " + str(nbSpecies)
        print(clusterCommand)
        ret = os.system(clusterCommand)
        if ret != 0:
            print("Command failed: ", ret)



        #clusterToNodenameCommand = "python3 " + scriptDir + "/unitigCluster_to_nodeName.py " + nodeToContigFilename + " "  ~/workspace/run/distance_matrix_composition.csv.gz.csv"

        if ret == 0: #Composition clustering can bug if distance matrix is empty (there is no contig of length >75k for isntance)
            clusterFilename = nodeToContigFilename.replace("Color", "Cluster")
            outputFile = open(clusterFilename, "w")
            outputFile.write("Name,Color\n")

            contigClusterFile = open(contigClusterFilename)
            contigClusterFile.readline()
            
            for line in contigClusterFile:
                line = line.rstrip()
                fields = line.split(",")

                unitigName = int(fields[0])
                cluster = int(fields[1])

                for nodeName in unitigToNodenames[unitigName]:
                    outputFile.write(str(nodeName) + "," + str(cluster) + "\n")

            outputFile.close()

    """
    
    return marker_contig_counts

def computeCompletenessContamination(marker_contig_counts):

    nbScgs = 0
    nbContaminatedScgs = 0

    for markerName, markerCounts in marker_contig_counts.items():
        nbScgs += 1
        if markerCounts[0] > 1: nbContaminatedScgs += 1

    completeness = float(nbScgs) / float(107)
    contamination = float(nbContaminatedScgs) / float(107)

    return completeness, contamination

def computeCompletenessContamination2(marker_contig_counts):

    nbScgs = 0
    nbContaminatedScgs = 0

    for markerName, markerCounts in marker_contig_counts.items():
        nbScgs += 1
        if markerCounts > 1: nbContaminatedScgs += 1

    completeness = float(nbScgs) / float(107)
    contamination = float(nbContaminatedScgs) / float(107)

    return completeness, contamination

    """
    logger.info("Number of contigs containing single-copy marker genes: " +
                str(len(contig_markers)))

    # Check if there are contigs with single-copy marker genes
    if len(contig_markers) == 0:
        logger.info(
            "Could not find contigs that contain single-copy marker genes. The dataset cannot be binned.")
        logger.info("Exiting MetaCoAG... Bye...!")
        sys.exit(1)

    # Get single-copy marker gene counts to make bins
    # -----------------------------------------------------

    logger.info("Determining contig counts for each single-copy marker gene")

    # Get the count of contigs for each single-copy marker gene
    my_gene_counts = list(marker_contig_counts.values())

    # Sort the counts in the descending order
    my_gene_counts.sort(reverse=True)

    logger.debug("Contig counts of single-copy marker genes: ")
    logger.debug(str(my_gene_counts))

    # Get contigs containing each single-copy marker gene for each iteration
    # -----------------------------------------------------

    smg_iteration = {}

    n = 0

    unique_my_gene_counts = list(set(my_gene_counts))
    unique_my_gene_counts.sort(reverse=True)

    # Get contigs for each iteration of single-copy marker gene
    for g_count in unique_my_gene_counts:

        # Get the single-copy marker genes with maximum count of contigs and
        # sort them in the descending order of the total marker genes contained

        total_contig_mgs = {}

        for item in marker_contig_counts:

            if marker_contig_counts[item] == g_count:

                total_contig_lengths = 0

                for contig in marker_contigs[item]:
                    contig_mg_counts = len(contig_markers[contig])
                    total_contig_lengths += contig_mg_counts

                total_contig_mgs[item] = total_contig_lengths

        total_contig_mgs_sorted = sorted(
            total_contig_mgs.items(), key=operator.itemgetter(1), reverse=True)

        for item in total_contig_mgs_sorted:
            smg_iteration[n] = marker_contigs[item[0]]
            n += 1

    # Initialise bins
    # -----------------------------------------------------

    bins = {}
    bin_of_contig = {}
    bin_markers = {}

    binned_contigs_with_markers = []

    logger.info("Initialising bins")

    # Initialise bins with the contigs having the first single-copy
    # marker gene according to the ordering

    for i in range(len(smg_iteration[0])):

        binned_contigs_with_markers.append(smg_iteration[0][i])
        contig_num = smg_iteration[0][i]

        bins[i] = [contig_num]
        bin_of_contig[contig_num] = i

        bin_markers[i] = contig_markers[contig_num]

    logger.debug("Number of initial bins detected: " +
                str(len(smg_iteration[0])))
    logger.debug("Initialised bins: ")
    logger.debug(bins)
    """

def computeScgClusters(contigFilename, nodeToContigFilename, mg_length_threshold, min_length, nbCores, contigName_to_contigIndex, contigLengths):

    unitigToNodenames = getUnitigToNodenames(nodeToContigFilename)




    scgSequencesFilename = contigFilename + ".scgSequences.ffn"
    #unitigSCG_filename = os.path.join(out_dir, genome_name + "_unitigSCG.fna")
    mmseqs_result_filename = contigFilename + ".scgSequences.mmseqs" #os.path.join(out_dir, genome_name + "_unitigSCG_mmseqs")
    mmseqs_tmp_filename = contigFilename + ".scgSequences.mmseqs.tmp" #os.path.join(out_dir, genome_name + "_unitigSCG_tmp")

    
    #print(unitigSCG_filename)
    command = "mmseqs easy-cluster " + scgSequencesFilename + " " + mmseqs_result_filename + " " + mmseqs_tmp_filename+ " --min-seq-id " + str(0.95) + " -c 0.8 --cov-mode 1 --alignment-mode 3 --threads 4"
    print(command)
    ret = os.system(command + " > /dev/null 2>&1")
    if ret != 0:
        print("MMseqs crash")
        sys.exit(ret)
    

    cluster_per_scg = {}
    cluster_per_scg, unitig_cluster_per_cog = count_cluster_per_SCG(contigFilename, cluster_per_scg, 0)

    scgIndex = 0
    scgName_to_index = {}
    for scgName in cluster_per_scg.keys():
        scgName_to_index[scgName] = scgIndex
        scgIndex += 1

    print(cluster_per_scg)
    print(unitig_cluster_per_cog)

    nodenamesUsed = {}
    #scgClusterName_to_index = {}
    #scgClusterIndex = 0


    scgClusterFile = open(nodeToContigFilename.replace("Color", "SCGcluster"), "w")
    scgClusterFile.write("Name,Color\n")
    scgClusterFile_contig = open(contigFilename.replace(".part", "") + ".scgCluster.txt", "w")
    #scgClusterFile.write("Name,Color\n")

    for scgName, scgClusters in unitig_cluster_per_cog.items():
        scgClusterIndexLocal = 0
        for scgCluster in scgClusters:
            scgClusterName = scgName + "_" + str(scgClusterIndexLocal)

            #if not scgClusterName in scgClusterName_to_index:
            #    scgClusterName_to_index[scgClusterName] = scgClusterIndex
            #    scgClusterIndex += 1

            for contigName in scgCluster:
                #print(contigName, scgClusterName)

                scgClusterFile_contig.write(contigName.replace("ctg", "") + "\t" + str(scgName_to_index[scgName]) + "\t" + str(scgClusterIndexLocal) + "\n")
                
                for nodeName in unitigToNodenames[contigName_to_contigIndex[contigName]]:
                    if nodeName in nodenamesUsed: continue

                    #scgClusterFile.write(str(nodeName) + "," + str(scgClusterName_to_index[scgClusterName]) + "\n")
                    scgClusterFile.write(str(nodeName) + "," + str(scgName_to_index[scgName]) + "_" + str(scgClusterIndexLocal) + "\n")
                    nodenamesUsed[nodeName] = True
                    break

            scgClusterIndexLocal += 1

    scgClusterFile.close()
    scgClusterFile_contig.close()

    return cluster_per_scg
    """
    cog_cluster_filename = contigFilename + ".testtmp"

    nbCogs = len(unitig_cluster_per_cog.keys())
    cog_index = 0
    unitig_to_cog = {}
    for cog_name in unitig_cluster_per_cog.keys():
        cluster_id = 1
        for unitig_cluster in unitig_cluster_per_cog[cog_name]:
            for unitig_name in unitig_cluster:

                if not unitig_name in unitig_to_cog:
                    unitig_to_cog[unitig_name] = [0] * nbCogs

                unitig_to_cog[unitig_name][cog_index] = cluster_id
                #print(unitig_name)
            cluster_id += 1
        cog_index += 1

    print(unitig_to_cog)

    cog_cluster_file = open(cog_cluster_filename, "w")
    header = "Name"
    for cog_name in unitig_cluster_per_cog.keys():
        header += "," + cog_name
    header += "\n"
    cog_cluster_file.write(header)

    for unitig_name in unitig_to_cog.keys():
        line = unitig_name
        for cluster in unitig_to_cog[unitig_name]:
            line += "," + str(cluster)
        line += "\n"
        cog_cluster_file.write(line)

    cog_cluster_file.close()
    """
    """
    scgClusterFile = open(nodeToContigFilename.replace("Color", "SCGcluster"), "w")
    scgClusterFile.write("Name,Color\n")

    nodenamesUsed = {}
    for unitigName, scgNames in contig_markers.items():
        for scgName in scgNames:
            for nodeName in unitigToNodenames[unitigName]:
                if nodeName in nodenamesUsed: continue

                scgClusterFile.write(str(nodeName) + "," + str(scgName_to_index[scgName]) + "\n")
                nodenamesUsed[nodeName] = True
                break

    scgClusterFile.close()
    """
    
def append_cluster(unitig_cluster, unitig_cluster_per_cog, cluster_per_scg, iteration):
    #print("-- break --")

    #print(cog_name)
    #if cog_name == "COG0081": lala += 1

    if unitig_cluster is not None:
        cog_name = unitig_cluster[0]
        if not cog_name in unitig_cluster_per_cog:
            unitig_cluster_per_cog[cog_name] = []
        unitig_cluster_per_cog[cog_name].append(unitig_cluster[1])
        #print("append cluster: ", cog_name, unitig_cluster[1])
        #if cog_name == "COG0541": print("append cluster: ", unitig_cluster[1])

        if not cog_name in cluster_per_scg:
            cluster_per_scg[cog_name] = []
        if iteration >= len(cluster_per_scg[cog_name]):
            cluster_per_scg[cog_name].append(0)
        cluster_per_scg[cog_name][iteration] += 1


def count_cluster_per_SCG(contigFilename, cluster_per_scg, iteration):

    nbUnitigs_per_cog = {}

    unitig_cluster_per_cog = {}
    unitig_cluster = None
    
    unitigSCG_cluster_filename = contigFilename + ".scgSequences.mmseqs_cluster.tsv" #os.path.join(out_dir, genome_name + "_unitigSCG_mmseqs_cluster.tsv")


    current_name = ""
    nbLines = 0
    #cluster_per_scg = {}
    i = 0
    for line in open(unitigSCG_cluster_filename):

        nbLines += 1
        line = line.rstrip()

        #print(line)
        name1, name2 = line.split("\t")
        unitig_name = name2.split("_")[0]
        
        
        cog_name = name1.split("_")[-1]

        #print(unitig_name, cog_name)
        #if cog_name == "COG0541":
        #    print(i, name1, name2, len(unitig_cluster), unitig_name)
        #    if "COG0541" in unitig_cluster_per_cog: print(unitig_cluster_per_cog["COG0541"])
        #    print(lili(unitig_cluster_per_cog))
        #    i +=1


        
        cog_name = name1.split("_")[-1]
        if not cog_name in nbUnitigs_per_cog:
            nbUnitigs_per_cog[cog_name] = {}
        nbUnitigs_per_cog[cog_name][name2.split("_")[0]] = 0
        #nbUnitigs_per_cog[cog_name][name1.split("_")[0]] = 0


        #loul=0
        #for cog_name in sorted(nbUnitigs_per_cog.keys()):
        #    loul += len(nbUnitigs_per_cog[cog_name])
        #print(line)
        #print(loul)
        #print(name1, name2)

        if name1 != current_name:
            append_cluster(unitig_cluster, unitig_cluster_per_cog, cluster_per_scg, iteration)

            cog_name = name1.split("_")[-1]
            unitig_cluster = (cog_name, [])

            current_name = name1

        unitig_cluster[1].append(unitig_name)


    append_cluster(unitig_cluster, unitig_cluster_per_cog, cluster_per_scg, iteration)
    #unitig_cluster_per_cog[cog_name].append(unitig_cluster)
    """
    for cog_name in sorted(unitig_cluster_per_cog.keys()):
        nbUnitigs = {}
        for cluster in unitig_cluster_per_cog[cog_name]:
            for unitig_name in cluster:
                nbUnitigs[unitig_name] = 0
        print(cog_name, len(nbUnitigs))
    """

    #print(cog_name, len(unitig_cluster_per_cog[cog_name]))

    #print(nbLines)
    #loul=0
    #for cog_name in sorted(nbUnitigs_per_cog.keys()):
    #    print(cog_name, len(nbUnitigs_per_cog[cog_name]))
    #    loul += len(nbUnitigs_per_cog[cog_name])
    #print(loul)

    return cluster_per_scg, unitig_cluster_per_cog





if __name__ == "__main__":
    main(sys.argv[1:])  
