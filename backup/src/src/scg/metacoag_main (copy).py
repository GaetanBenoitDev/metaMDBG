#!/usr/bin/env python3

import sys
import time
import argparse
import os
import math
import logging
import operator
import gc
import subprocess
import pathlib
import concurrent.futures
import csv

from Bio import SeqIO
from igraph import *
from networkx.classes.function import set_node_attributes
from tqdm import tqdm

from metacoag_utils import feature_utils
from metacoag_utils import marker_gene_utils
from metacoag_utils import matching_utils
from metacoag_utils import label_prop_utils
from metacoag_utils import graph_utils
from metacoag_utils.bidirectionalmap import BidirectionalMap

# Set paramters
# ---------------------------------------------------

MAX_WEIGHT = sys.float_info.max
M_MARKER_GENES = 108


# Setup argument parser
# ---------------------------------------------------

ap = argparse.ArgumentParser(description="""MetaCoAG is a NGS data-based metagenomic contig binning tool that makes use of the 
connectivity information found in assembly graphs, apart from the composition and coverage information. 
MetaCoAG makes use of single-copy marker genes along with a graph matching technique and a label propagation technique to bin contigs.""")

ap.add_argument("--assembler", required=True,
                help="name of the assembler used")
ap.add_argument("--contigs", required=True, help="path to the contigs file")
ap.add_argument("--graph", required=True,
                help="path to the assembly graph file")
ap.add_argument("--paths", required=True,
                help="path to the contigs.paths file")
ap.add_argument("--abundance", required=True,
                help="path to the abundance file")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, default='',
                help="prefix for the output file")
ap.add_argument("--min_length", required=False, type=int, default=1000,
                help="minimum length of contigs to consider for compositional probability. [default: 1000]")
ap.add_argument("--p_intra", required=False, type=float, default=0.1,
                help="minimum probability of an edge matching to assign to the same bin. [default: 0.1]")
ap.add_argument("--p_inter", required=False, type=float, default=0.01,
                help="maximum probability of an edge matching to create a new bin. [default: 0.01]")
ap.add_argument("--depth", required=False, type=int, default=10,
                help="depth to consider for label propagation. [default: 10]")
ap.add_argument("--mg_threshold", required=False, type=float, default=0.5,
                help="length threshold to consider marker genes. [default: 0.5]")
ap.add_argument("--bin_mg_threshold", required=False, type=float, default=0.33333,
                help="minimum fraction of marker genes that should be present in a bin. [default: 0.33333]")
ap.add_argument("--min_bin_size", required=False, type=int, default=200000,
                help="minimum size of a bin to output in base pairs. [default: 200000]")
ap.add_argument("--hmm", required=False, type=str, default="",
                help="path to marker.hmm file. [default: auxiliary/marker.hmm]")
ap.add_argument("--d_limit", required=False, type=int, default=20,
                help="distance limit for contig matching. [default: 20]")
ap.add_argument("--delimiter", required=False, type=str, default=",",
                help="delimiter for output results. [default: , (comma)]")
ap.add_argument("--nthreads", required=False, type=int, default=8,
                help="number of threads to use. [default: 8]")

# Parse arguments
args = vars(ap.parse_args())
assembler = args["assembler"]
contigs_file = args["contigs"]
assembly_graph_file = args["graph"]
contig_paths_file = args["paths"]
abundance_file = args["abundance"]
output_path = args["output"]
prefix = args["prefix"]
min_length = args["min_length"]
p_intra = args["p_intra"]
p_inter = args["p_inter"]
depth = args["depth"]
mg_threshold = args["mg_threshold"]
bin_mg_threshold = args["bin_mg_threshold"]
min_bin_size = args["min_bin_size"]
hmm = args["hmm"]
d_limit = args["d_limit"]
delimiter = args["delimiter"]
nthreads = args["nthreads"]

if hmm == "":
    hmm = os.path.join(pathlib.Path(
        __file__).parent.absolute().parent, 'auxiliary', 'marker.hmm')

bin_threshold = -math.log(p_intra, 10)
break_threshold = -math.log(p_inter, 10)

n_bins = 0


# Setup logger
# -----------------------

logger = logging.getLogger('MetaCoaAG 1.0')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
consoleHeader.setLevel(logging.INFO)
logger.addHandler(consoleHeader)

# Setup output path for log file
fileHandler = logging.FileHandler(output_path + "/" + prefix + "metacoag.log")
fileHandler.setLevel(logging.DEBUG)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

logger.info(
    "Welcome to MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs.")
logger.info("This version of MetaCoAG makes use of the assembly graph produced by SPAdes which is based on the de Bruijn graph approach.")

logger.info("Input arguments: ")
logger.info("Assembler used: " + assembler)
logger.info("Contigs file: " + contigs_file)
logger.info("Assembly graph file: " + assembly_graph_file)
logger.info("Contig paths file: " + contig_paths_file)
logger.info("Abundance file: " + abundance_file)
logger.info("Final binning output file: " + output_path)
logger.info("Marker file: " + hmm)
logger.info("Minimum length of contigs to consider: " + str(min_length))
logger.info("Depth to consider for label propagation: " + str(depth))
logger.info("p_intra: " + str(p_intra))
logger.info("p_inter: " + str(p_inter))
logger.debug("bin_threshold: " + str(bin_threshold))
logger.debug("break_threshold: " + str(break_threshold))
logger.info("mg_threshold: " + str(mg_threshold))
logger.info("bin_mg_threshold: " + str(bin_mg_threshold))
logger.info("min_bin_size: " + str(min_bin_size) + " base pairs")
logger.info("d_limit: " + str(d_limit))
logger.info("Number of threads: " + str(nthreads))

logger.info("MetaCoAG started")

start_time = time.time()


# Get links of the assembly graph
# ------------------------------------------------------------------------

try:
    if assembler == "spades":

        # Get paths, segments, links and contigs of the assembly graph
        paths, segment_contigs, contig_segments, node_count, contigs_map, contig_names = graph_utils.get_segment_paths_spades(
            contig_paths_file)

        # Get reverse mapping of contig map
        contigs_map_rev = contigs_map.inverse

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

    if assembler == "megahit":

        original_contigs = {}
        contig_descriptions = {}

        # Get mapping of original contig identifiers with descriptions
        for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
            original_contigs[record.id] = str(record.seq)
            contig_descriptions[record.id] = record.description

        # Get links and contigs of the assembly graph
        node_count, graph_contigs, links, contig_names = graph_utils.get_links_megahit(
            assembly_graph_file)

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

    if assembler == "flye":

        # Get links and contigs of the assembly graph
        node_count, links, contig_names = graph_utils.get_links_flye(
            assembly_graph_file)

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

except:
    logger.error(
        "Please make sure that the correct path to the contig paths file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")
    sys.exit(1)


# Construct the assembly graph
# -------------------------------

all_contigs = [x for x in range(node_count)]

try:

    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)
    logger.info("Total number of contigs available: " + str(node_count))

    # Name vertices with contig identifiers
    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["label"] = contig_names[i]

    # Get list of edges
    if assembler == "spades":
        edge_list = graph_utils.get_graph_edges_spades(
            assembly_graph_file=assembly_graph_file,
            contigs_map=contigs_map,
            contigs_map_rev=contigs_map_rev,
            paths=paths,
            segment_contigs=segment_contigs)

    if assembler == "flye":
        edge_list = graph_utils.get_graph_edges_flye(
            links=links,
            contig_names_rev=contig_names_rev)

    if assembler == "megahit":
        edge_list = graph_utils.get_graph_edges_megahit(
            links=links,
            contig_names_rev=contig_names_rev)

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)

    # Simplify the graph
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    logger.info("Total number of edges in the assembly graph: " +
                str(len(list(assembly_graph.es))))

except:
    logger.error(
        "Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting MetaCoAG... Bye...!")
    sys.exit(1)


if assembler == "megahit":

    # Map original contig identifiers to contig identifiers of MEGAHIT assembly graph
    graph_to_contig_map = BidirectionalMap()

    for (n, m), (n2, m2) in zip(graph_contigs.items(), original_contigs.items()):
        if m == m2:
            graph_to_contig_map[n] = n2

    graph_to_contig_map_rev = graph_to_contig_map.inverse


# Get the number of samples and the length and coverage of contigs
# ------------------------------------------------------------------------

logger.info("Obtaining lengths and coverage values of contigs")

if assembler == "megahit":
    sequences, coverages, contig_lengths, n_samples = feature_utils.get_cov_len_megahit(
        contigs_file=contigs_file,
        contig_names_rev=contig_names_rev,
        graph_to_contig_map_rev=graph_to_contig_map_rev,
        min_length=min_length,
        abundance_file=abundance_file)

else:
    sequences, coverages, contig_lengths, n_samples = feature_utils.get_cov_len(
        contigs_file=contigs_file,
        contig_names_rev=contig_names_rev,
        min_length=min_length,
        abundance_file=abundance_file)


# Set intra weight and inter weight
# ------------------------------------------------------------------------

w_intra = bin_threshold * (n_samples + 1)
w_inter = break_threshold * (n_samples + 1)

logger.debug("w_intra: " + str(w_intra))
logger.debug("w_inter: " + str(w_inter))


# Get tetramer composition of contigs
# ------------------------------------------------------------------------

logger.info("Obtaining tetranucleotide frequencies of contigs")

normalized_tetramer_profiles = feature_utils.get_tetramer_profiles(
    output_path=output_path,
    sequences=sequences,
    contig_lengths=contig_lengths,
    min_length=min_length,
    nthreads=nthreads)

del sequences
gc.collect()


# Get contigs with marker genes
# -----------------------------------------------------

logger.info("Scanning for single-copy marker genes")

if not os.path.exists(contigs_file + ".hmmout"):
    # Run FragGeneScan and HMMER if .hmmout file is not present
    logger.info("Obtaining hmmout file")
    marker_gene_utils.scan_for_marker_genes(
        contigs_file=contigs_file,
        nthreads=nthreads,
        markerURL=hmm)
else:
    logger.info(".hmmout file already exists")

logger.info("Obtaining contigs with single-copy marker genes")

# Get contigs with single-copy marker genes and count of contigs for each single-copy marker gene
if assembler == "megahit":
    marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes_megahit(
        contigs_file=contigs_file,
        contig_names_rev=contig_names_rev,
        graph_to_contig_map_rev=graph_to_contig_map_rev,
        mg_length_threshold=mg_threshold,
        contig_lengths=contig_lengths,
        min_length=min_length)

else:
    marker_contigs, marker_contig_counts, contig_markers = marker_gene_utils.get_contigs_with_marker_genes(
        contigs_file=contigs_file,
        contig_names_rev=contig_names_rev,
        mg_length_threshold=mg_threshold,
        contig_lengths=contig_lengths,
        min_length=min_length)

    all_contig_markers = marker_gene_utils.get_all_contigs_with_marker_genes(
        contigs_file=contigs_file,
        contig_names_rev=contig_names_rev,
        mg_length_threshold=mg_threshold)

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

# Assign contigs with single-copy marker genes to bins
# -----------------------------------------------------

logger.info(
    "Matching and assigning contigs with single-copy marker genes to bins")

# Matching contigs to bins
bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers = matching_utils.match_contigs(
    smg_iteration=smg_iteration,
    bins=bins,
    n_bins=n_bins,
    bin_of_contig=bin_of_contig,
    binned_contigs_with_markers=binned_contigs_with_markers,
    bin_markers=bin_markers,
    contig_markers=contig_markers,
    contig_lengths=contig_lengths,
    contig_names=contig_names,
    normalized_tetramer_profiles=normalized_tetramer_profiles,
    coverages=coverages,
    assembly_graph=assembly_graph,
    w_intra=w_intra,
    w_inter=w_inter,
    d_limit=d_limit)

logger.debug("Number of bins after matching: " + str(len(bins)))

logger.debug("Bins with contigs containing seed marker genes")

for b in bins:
    logger.debug(str(b) + ": " + str(bins[b]))

logger.debug(
    "Number of binned contigs with single-copy marker genes: " + str(len(bin_of_contig)))

del smg_iteration
del my_gene_counts
del unique_my_gene_counts
del marker_contigs
del marker_contig_counts
del total_contig_mgs
gc.collect()

# Further assign contigs with seed marker genes
# -------------------------------------------------

# Get contigs with single-copy marker genes which are not matched to bins
unbinned_mg_contigs = list(
    set(contig_markers.keys()) - set(binned_contigs_with_markers))

unbinned_mg_contig_lengths = {}

# Get the lengths of the unmatched contigs
for contig in unbinned_mg_contigs:
    contigid = contig
    unbinned_mg_contig_lengths[contig] = contig_lengths[contigid]

# Sort the unmatched in the the descending order of their contig lengths
unbinned_mg_contig_lengths_sorted = sorted(
    unbinned_mg_contig_lengths.items(), key=operator.itemgetter(1), reverse=True)

logger.debug("Number of unbinned contigs with single-copy marker genes: " +
             str(len(unbinned_mg_contigs)))

logger.info("Further assigning contigs with single-copy marker genes")

# Further assigning unmatched contigs to bins
bins, bin_of_contig, n_bins, bin_markers, binned_contigs_with_markers = matching_utils.further_match_contigs(
    unbinned_mg_contigs=unbinned_mg_contig_lengths_sorted,
    min_length=min_length,
    bins=bins,
    n_bins=n_bins,
    bin_of_contig=bin_of_contig,
    binned_contigs_with_markers=binned_contigs_with_markers,
    bin_markers=bin_markers,
    contig_markers=contig_markers,
    normalized_tetramer_profiles=normalized_tetramer_profiles,
    coverages=coverages,
    w_intra=w_intra)

# Get remaining contigs with single-copy marker genes which are not assigned to bins
unbinned_mg_contigs = list(
    set(contig_markers.keys()) - set(binned_contigs_with_markers))

logger.debug("Remaining number of unbinned MG seed contigs: " +
             str(len(unbinned_mg_contigs)))
logger.debug(
    "Number of binned contigs with single-copy marker genes: " + str(len(bin_of_contig)))

del unbinned_mg_contigs
del unbinned_mg_contig_lengths
del unbinned_mg_contig_lengths_sorted
gc.collect()

# Get seed bin counts and profiles
# -----------------------------------------------------

smg_bin_counts = []

# Get the count of contigs with single-copy marker genes in each bin
for i in bins:
    smg_bin_counts.append(len(bins[i]))

# Get composition and coverage profiles for each bin
# based on contigs with single-copy marker genes
bin_seed_tetramer_profiles, bin_seed_coverage_profiles = feature_utils.get_bin_profiles(
    bins=bins,
    coverages=coverages,
    normalized_tetramer_profiles=normalized_tetramer_profiles)

# Get binned and unbinned contigs
# -----------------------------------------------------

binned_contigs = list(bin_of_contig.keys())

unbinned_contigs = list(
    set([x for x in range(node_count)]) - set(binned_contigs))

logger.debug("Number of binned contigs: " + str(len(binned_contigs)))
logger.debug("Number of unbinned contigs: " + str(len(unbinned_contigs)))
logger.debug("Number of binned contigs with markers: " +
             str(len(binned_contigs_with_markers)))

# Get isolated vertices and components without labels
# -----------------------------------------------------

# Get isolated contigs with no neighbours
isolated = graph_utils.get_isolated(node_count, assembly_graph)

# Get connected contigs within the labelled components
non_isolated = graph_utils.get_non_isolated(
    node_count=node_count,
    assembly_graph=assembly_graph,
    binned_contigs=binned_contigs)

logger.debug("Number of non-isolated contigs: " + str(len(non_isolated)))

# Propagate labels to vertices of unlabelled long contigs
# -----------------------------------------------------

logger.info(
    "Propagating labels to connected vertices of unlabelled long contigs")

# Label propagation on connected vertices of unlabelled long contigs
bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.label_prop(
    bin_of_contig=bin_of_contig,
    bins=bins,
    contig_markers=contig_markers,
    bin_markers=bin_markers,
    binned_contigs_with_markers=binned_contigs_with_markers,
    smg_bin_counts=smg_bin_counts,
    non_isolated=non_isolated,
    contig_lengths=contig_lengths,
    min_length=min_length,
    assembly_graph=assembly_graph,
    normalized_tetramer_profiles=normalized_tetramer_profiles,
    coverages=coverages,
    depth=1,
    weight=w_intra)

logger.debug("Total number of binned contigs: " + str(len(bin_of_contig)))


# Further propagate labels to vertices of unlabelled long contigs
# --------------------------------------------------------------------------------

# Further label propagation on connected vertices of unlabelled long contigs
bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.label_prop(
    bin_of_contig=bin_of_contig,
    bins=bins,
    contig_markers=contig_markers,
    bin_markers=bin_markers,
    binned_contigs_with_markers=binned_contigs_with_markers,
    smg_bin_counts=smg_bin_counts,
    non_isolated=non_isolated,
    contig_lengths=contig_lengths,
    min_length=min_length,
    assembly_graph=assembly_graph,
    normalized_tetramer_profiles=normalized_tetramer_profiles,
    coverages=coverages,
    depth=depth,
    weight=w_inter)

logger.debug("Total number of binned contigs: " + str(len(bin_of_contig)))


# Get binned and unbinned contigs
# -----------------------------------------------------

binned_contigs = list(bin_of_contig.keys())

unbinned_contigs = list(
    set([x for x in range(node_count)]) - set(binned_contigs))

logger.debug("Number of binned contigs: " + str(len(binned_contigs)))
logger.debug("Number of unbinned contigs: " + str(len(unbinned_contigs)))


# Propagate labels to vertices of unlabelled long contigs in isolated components
# -----------------------------------------------------------------------------------------------

logger.info(
    "Further propagating labels to vertices of unlabelled long contigs")

# Get long unbinned contigs
long_unbinned = list(
    filter(lambda contig: contig not in bin_of_contig and contig_lengths[contig] >= min_length, all_contigs))

# Starting propagation of labels to vertices of unlabelled long contigs

assigned = [None for itr in long_unbinned]

executor = concurrent.futures.ThreadPoolExecutor(max_workers=nthreads)

# Thread function for workers


def thread_function(n, contig, coverages, normalized_tetramer_profiles, bin_seed_tetramer_profiles, bin_seed_coverage_profiles):
    bin_result = label_prop_utils.assign_long(
        contigid=contig,
        coverages=coverages,
        normalized_tetramer_profiles=normalized_tetramer_profiles,
        bin_tetramer_profiles=bin_seed_tetramer_profiles,
        bin_coverage_profiles=bin_seed_coverage_profiles)
    assigned[n] = bin_result


exec_args = []

for n, contig in enumerate(long_unbinned):
    exec_args.append(
        (n, contig, coverages, normalized_tetramer_profiles, bin_seed_tetramer_profiles, bin_seed_coverage_profiles))

for itr in tqdm(executor.map(lambda p: thread_function(*p), exec_args), total=len(long_unbinned)):
    pass

executor.shutdown(wait=True)

# End propagation of labels to vertices of unlabelled long contigs

put_to_bins = [x for x in assigned if x is not None]

if len(put_to_bins) == 0:
    logger.debug("No further contigs were binned")
else:
    # Add contigs to bins according to assignment
    bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.assign_to_bins(
        put_to_bins=put_to_bins,
        bins=bins,
        bin_of_contig=bin_of_contig,
        bin_markers=bin_markers,
        binned_contigs_with_markers=binned_contigs_with_markers,
        contig_markers=contig_markers,
        contig_lengths=contig_lengths)

logger.debug("Total number of binned contigs: " + str(len(bin_of_contig)))

#  Further propagate labels to vertices of unlabelled long contigs
# --------------------------------------------------------------------------------

logger.info(
    "Further propagating labels to connected vertices of unlabelled long contigs")

# Further label propagation on connected vertices of unlabelled long contigs
bins, bin_of_contig, bin_markers, binned_contigs_with_markers = label_prop_utils.final_label_prop(
    bin_of_contig=bin_of_contig,
    bins=bins,
    contig_markers=contig_markers,
    bin_markers=bin_markers,
    binned_contigs_with_markers=binned_contigs_with_markers,
    smg_bin_counts=smg_bin_counts,
    contig_lengths=contig_lengths,
    min_length=min_length,
    assembly_graph=assembly_graph,
    normalized_tetramer_profiles=normalized_tetramer_profiles,
    coverages=coverages,
    depth=depth,
    weight=MAX_WEIGHT)

logger.debug("Total number of binned contigs: " + str(len(bin_of_contig)))


# Get elapsed time
# -----------------------------------

# Determine elapsed time
elapsed_time = time.time() - start_time

# Print elapsed time for the process
logger.info("Elapsed time: " + str(elapsed_time) + " seconds")


# Get bin sizes
# -----------------------------------

bin_size = {}

for b in bins:
    bin_size[b] = 0
    for contig in bins[b]:
        bin_size[b] += contig_lengths[contig]


# Merge bins
# -----------------------------------

# Create graph
bins_graph = Graph()

# Add vertices
bins_graph.add_vertices(len(bins))

# Name vertices with contig identifiers
for i in range(len(bins)):
    bins_graph.vs[i]["id"] = i
    bins_graph.vs[i]["label"] = "bin " + str(i+1)

bins_to_rem = []

for b in bins:

    possible_bins = []

    no_possible_bins = True

    logger.debug("Bin " + str(b) + ": # contigs: " + str(len(bins[b])) + ", bin size: " + str(
        bin_size[b]) + "bp, # markers: " + str(len(bin_markers[b])))

    min_pb = -1
    min_pb_weight = MAX_WEIGHT

    for pb in bin_markers:
        common_mgs = list(set(bin_markers[pb]).intersection(
            set(bin_markers[b])))

        if len(common_mgs) == 0:

            tetramer_dist = matching_utils.get_tetramer_distance(bin_seed_tetramer_profiles[b],
                                                                 bin_seed_tetramer_profiles[pb])
            prob_comp = matching_utils.get_comp_probability(
                tetramer_dist)
            prob_cov = matching_utils.get_cov_probability(
                bin_seed_coverage_profiles[pb], bin_seed_coverage_profiles[b])
            prob_product = prob_comp * prob_cov
            log_prob = 0

            if prob_product > 0.0:
                log_prob = - (math.log(prob_comp, 10) +
                              math.log(prob_cov, 10))
            else:
                log_prob = MAX_WEIGHT

            if log_prob <= w_intra:

                prob_cov1 = matching_utils.get_cov_probability(
                    bin_seed_coverage_profiles[pb], bin_seed_coverage_profiles[b])
                prob_product1 = prob_comp * prob_cov1
                log_prob1 = 0

                if prob_product1 > 0.0:
                    log_prob1 = - \
                        (math.log(prob_comp, 10) + math.log(prob_cov1, 10))
                else:
                    log_prob1 = MAX_WEIGHT

                if log_prob1 <= w_intra:

                    possible_bins.append(pb)

                    if log_prob < min_pb_weight:
                        min_pb_weight = log_prob
                        min_pb = pb

    if min_pb != -1:
        bins_graph.add_edge(b, min_pb)
        no_possible_bins = False

    if no_possible_bins and len(bin_markers[b]) < M_MARKER_GENES * bin_mg_threshold:
        bins_to_rem.append(b)

bin_cliques = bins_graph.maximal_cliques()


# Get bin clique sizes
# -----------------------------------

bin_clique_size = {}

for bin_clique in bin_cliques:

    bin_name = '_'.join(str(x) for x in list(bin_clique))

    bin_clique_size[bin_name] = 0

    for b in bin_clique:
        bin_clique_size[bin_name] += bin_size[b]


# Get final list of bins and write result to output file
# ----------------------------------------------------------

# Get output path
output_bins_path = output_path + prefix + "bins/"
lq_output_bins_path = output_path + prefix + "low_quality_bins/"

# Create output directory for bin files
if not os.path.isdir(output_bins_path):
    subprocess.run("mkdir -p " + output_bins_path, shell=True)
if not os.path.isdir(lq_output_bins_path):
    subprocess.run("mkdir -p " + lq_output_bins_path, shell=True)

final_bins = {}
lowq_bins = {}

final_bin_count = 0

with open(output_path + prefix + "contig_to_bin.tsv", mode='w') as out_file:
    output_writer = csv.writer(
        out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)

    for bin_clique in bin_cliques:

        bin_name = '_'.join(str(x) for x in list(bin_clique))

        can_write = True

        if len(bin_clique) == 1 and bin_clique[0] in bins_to_rem:
            can_write = False

        if can_write and bin_clique_size[bin_name] >= min_bin_size:

            final_bin_count += 1

            for b in bin_clique:

                for contig in bins[b]:
                    final_bins[contig] = bin_name

                    if assembler == "megahit":
                        output_writer.writerow(
                            [contig_descriptions[graph_to_contig_map[contig_names[contig]]], "bin_" + bin_name])
                    else:
                        output_writer.writerow(
                            [contig_names[contig], "bin_" + bin_name])

        else:

            for b in bin_clique:

                for contig in bins[b]:
                    lowq_bins[contig] = bin_name


logger.info("Writing the Final Binning result to file")

bin_files = {}

for bin_name in set(final_bins.values()):
    bin_files[bin_name] = open(
        output_bins_path + prefix + "bin_" + bin_name + ".fasta", 'w+')

for bin_name in set(lowq_bins.values()):
    bin_files[bin_name] = open(
        lq_output_bins_path + prefix + "bin_" + bin_name + "_seqs.fasta", 'w+')


for n, record in tqdm(enumerate(SeqIO.parse(contigs_file, "fasta")), desc="Splitting contigs into bins"):

    if assembler == "megahit":
        contig_num = contig_names_rev[graph_to_contig_map_rev[record.id]]
    else:
        contig_num = contig_names_rev[record.id]

    if contig_num in final_bins:
        bin_files[final_bins[contig_num]].write(
            f'>{str(record.id)}\n{str(record.seq)}\n')

    elif contig_num in lowq_bins:
        bin_files[lowq_bins[contig_num]].write(
            f'>{str(record.id)}\n{str(record.seq)}\n')

# Close output files
for c in set(final_bins.values()):
    bin_files[c].close()

for c in set(lowq_bins.values()):
    bin_files[c].close()


logger.info("Producing " + str(final_bin_count) + " bins...")
logger.info("Final binning results can be found in " + str(output_bins_path))


# Exit program
# -----------------------------------

logger.info("Thank you for using MetaCoAG!")
