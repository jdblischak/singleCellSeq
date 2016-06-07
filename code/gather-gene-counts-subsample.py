#!/usr/bin/env python

# This file is a variant of gather-gene-counts.py used in the main
# analysis pipeline. Since the processing for the subsampling analysis
# is sufficiently different, I reasoned it would be easier (and less
# error prone) to just create another script instead of complicating
# the logic of the original to handle the disparate use cases.

# In the main analysis, the per sample reads do not have their own
# genecounts file. This is because they can be accurately summed from
# their per lane counts. Not creating the "combined" counts files for
# the reads reduces the amount of computation and storage
# requirements. In contrast, we do not need per lane counts for the
# subsampling analyses. Thus it is more efficient to only work with
# the per sample ("combined") counts.

help_message = """\
Gathers the gene counts for each sample for the subsampling analysis.

Usage:
    gather-gene-counts-subsample.py prefix [files]

    prefix - prefix for names of output files,
             e.g. outdir/ or outdir/250000-,
             outdir must already exist
    files - At least one featureCounts output filename,
            e.g. counts/*.250000.*genecounts.txt
    """

# Notes:
#   * Should be run from data directory
#   * This file can be run with either Python 2 or 3

# Creates 2 files:
#
# 1. reads-raw-single-per-sample.txt
# 2. molecules-raw-single-per-sample.txt

# Explanation of file names:
#
# reads vs. molecules:
#     - reads are number of sequences per gene, as per traditional RNA-seq
#     - molecules are the number of UMIs per gene

################################################################################

import sys

if "-h" in sys.argv or "--help" in sys.argv:
    sys.exit(help_message)
elif len(sys.argv) < 3:
    sys.exit("Error: Need at least two arguments.\n\n%s"%(help_message))
else:
    prefix = sys.argv[1]
    files = sys.argv[2:]

# For testing:
# import glob
# files = glob.glob("counts/*genecounts.txt")
# print(len(files))
# sys.exit()

################################################################################
# Collate subsampled per sample reads and per sample molecules for single cells
################################################################################

reads_raw_single_per_sample = open(prefix + "reads-raw-single-per-sample.txt", "w")
molecules_raw_single_per_sample = open(prefix + "molecules-raw-single-per-sample.txt", "w")

# Get gene names from first file
gene_list = []
f = files[0]
handle = open(f, "r")
for line in handle:
    if line[0] == "#" or line[:6] == "Geneid":
        continue
    cols = line.strip("\n").split("\t")
    gene = cols[0]
    gene_list.append(gene)
handle.close()

# Create headers
per_sample_header = "individual\treplicate\twell\t"
reads_raw_single_per_sample.write(per_sample_header + "\t".join(gene_list) + "\n")
molecules_raw_single_per_sample.write(per_sample_header + "\t".join(gene_list) + "\n")

# Process each input genecounts file
for f in files:

    # Get counts from f
    g = 0 # iterator for indexing gene names
    counts = ["NA"] * len(gene_list)
    handle = open(f, "r")
    for line in handle:
        if line[0] == "#" or line[:6] == "Geneid":
            continue
        cols = line.strip("\n").split("\t")
        gene = cols[0]
        assert gene == gene_list[g], "Gene names are not in correct order in file: %s"%(f)
        counts[g] = cols[6]
        g += 1
    handle.close()

    # Get meta data from filename
    dir, fname = f.split("/")
    fname_parts = fname.split(".")
    individual, replicate, well = fname_parts[:3]
    individual = "NA" + individual
    replicate = "r" + replicate

    # molecules-raw-single-per-sample.txt
    if "combined" in fname and "rmdup" in fname:
        molecules_raw_single_per_sample.write(
          individual + "\t" + replicate + "\t" + well + "\t" + \
          "\t".join(counts) + "\n")
    # reads-raw-single-per-sample.txt
    elif "combined" in fname and "rmdup" not in fname:
        reads_raw_single_per_sample.write(
          individual + "\t" + replicate + "\t" + well + "\t" + \
          "\t".join(counts) + "\n")
    else:
        sys.stderr.write("Incompatible filename. Skipped %s\n"%(f))

# Close connections to files
reads_raw_single_per_sample.close()
molecules_raw_single_per_sample.close()
