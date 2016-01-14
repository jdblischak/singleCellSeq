#!/usr/bin/env python

help_message = """\
Gathers the gene counts for each sample.

Usage:
    gather-gene-counts.py prefix [files]

    prefix - prefix for names of output files,
             e.g. outdir/ or outdir/name-
    files - At least one featureCounts output filename,
            e.g. counts/*genecounts.txt
    """

# Notes:
#   * Should be run from data directory
#   * This file can be run with either Python 2 or 3

# Creates 6 files:
#
# 1. reads-raw-bulk-per-lane.txt
# 2. reads-raw-bulk-per-sample.txt
# 3. reads-raw-single-per-lane.txt
# 4. reads-raw-single-per-sample.txt
# 5. molecules-raw-single-per-lane.txt
# 6. molecules-raw-single-per-sample.txt

# Explanation of file names:
#
# reads vs. molecules:
#     - reads are number of sequences per gene, as per traditional RNA-seq
#     - molecules are the number of UMIs per gene
# bulk vs. single:
#     - bulk is sequencing of a population of cells, as per traditional RNA-seq
#     - single is sequencing of single cells
# lane vs. sample:
#     - lane is the gene counts for a given sample from one lane of sequencing
#     - sample is the sum of the gene counts from all lanes for a given sample

################################################################################

import sys
import pandas as pd

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
# Collate per lane reads and per sample molecules
################################################################################

# Files 1, 3, 5, and 6 can be created directly from the data.

reads_raw_bulk_per_lane = open(prefix + "reads-raw-bulk-per-lane.txt", "w")
reads_raw_single_per_lane = open(prefix + "reads-raw-single-per-lane.txt", "w")
molecules_raw_single_per_lane = open(prefix + "molecules-raw-single-per-lane.txt", "w")
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
per_lane_header = per_sample_header + "index\tlane\tflow_cell\t"
reads_raw_bulk_per_lane.write(per_lane_header + "\t".join(gene_list) + "\n")
reads_raw_single_per_lane.write(per_lane_header + "\t".join(gene_list) + "\n")
molecules_raw_single_per_lane.write(per_lane_header + "\t".join(gene_list) + "\n")
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

    # 6. molecules-raw-single-per-sample.txt
    if "combined" in fname:
        molecules_raw_single_per_sample.write(
          individual + "\t" + replicate + "\t" + well + "\t" + \
          "\t".join(counts) + "\n")
    else:
        index, lane = fname_parts[3:5]
        flow_cell = fname_parts[6]
        # 5. molecules-raw-single-per-lane.txt
        if "rmdup" in fname:
            molecules_raw_single_per_lane.write(
              individual + "\t" + replicate + "\t" + well + "\t" + \
              index + "\t" + lane + "\t" + flow_cell + "\t" + \
              "\t".join(counts) + "\n")
        else:
            # 1. reads-raw-bulk-per-lane.txt
            if "bulk" in fname:
                reads_raw_bulk_per_lane.write(
                  individual + "\t" + replicate + "\t" + well + "\t" + \
                  index + "\t" + lane + "\t" + flow_cell + "\t" + \
                   "\t".join(counts) + "\n")
            else:
                # 3. reads-raw-single-per-lane.txt
                reads_raw_single_per_lane.write(
                  individual + "\t" + replicate + "\t" + well + "\t" + \
                  index + "\t" + lane + "\t" + flow_cell + "\t" + \
                  "\t".join(counts) + "\n")


# Close connections to files
reads_raw_bulk_per_lane.close()
reads_raw_single_per_lane.close()
molecules_raw_single_per_lane.close()
molecules_raw_single_per_sample.close()

################################################################################
# Sum reads per lane to obtain reads per sample
################################################################################

# 2. reads-raw-bulk-per-sample.txt
reads_raw_bulk_per_lane_df = pd.read_table(prefix + "reads-raw-bulk-per-lane.txt")
reads_raw_bulk_per_sample_df = \
    reads_raw_bulk_per_lane_df.groupby(["individual", "replicate", "well"],
    as_index = False).sum()
reads_raw_bulk_per_sample_df.to_csv(prefix + "reads-raw-bulk-per-sample.txt",
    sep = "\t", na_rep = "NA", index = False)

# 4. reads-raw-single-per-sample.txt
reads_raw_single_per_lane_df = pd.read_table(prefix + "reads-raw-single-per-lane.txt")
reads_raw_single_per_sample_df = \
    reads_raw_single_per_lane_df.groupby(["individual", "replicate", "well"],
    as_index = False).sum()
reads_raw_single_per_sample_df.to_csv(prefix + "reads-raw-single-per-sample.txt",
    sep = "\t", na_rep = "NA", index = False)

################################################################################
# Test expected number of rows of each file
################################################################################

def check_line_num(fname, num_lines):
    """
    Assert that the file fname has num_lines lines.
    """
    handle = open(fname, "r")
    line_count = 0
    for line in handle:
        line_count += 1
    assert line_count == num_lines, "%s has %d lines, not %d"%(
      fname, line_count, num_lines)

# 1. reads-raw-bulk-per-lane.txt
# *  bulk reads per lane: 3 individuals * 3 reps * 2 indexes * 4 lanes = 72
check_line_num(prefix + "reads-raw-bulk-per-lane.txt", 72 + 1)

# 2. reads-raw-bulk-per-sample.txt
# *  bulk reads per sample: 3 individuals * 3 reps = 9
check_line_num(prefix + "reads-raw-bulk-per-sample.txt", 9 + 1)

# 3. reads-raw-single-per-lane.txt
# *  single cell reads per lane: 3 individuals * 3 reps * 96 wells * 3 lanes = 2,592
check_line_num(prefix + "reads-raw-single-per-lane.txt", 2592 + 1)

# 4. reads-raw-single-per-sample.txt
# *  single cell reads per sample: 3 individuals * 3 reps * 96 wells = 864
check_line_num(prefix + "reads-raw-single-per-sample.txt", 864 + 1)

# 5. molecules-raw-single-per-lane.txt
# *  single cell molecules per lane: 3 individuals * 3 reps * 96 wells * 3 lanes = 2,592
check_line_num(prefix + "molecules-raw-single-per-lane.txt", 2592 + 1)

# 6. molecules-raw-single-per-sample.txt
# *  single cell molecules per sample: 3 individuals * 3 reps * 96 wells = 864
check_line_num(prefix + "molecules-raw-single-per-sample.txt", 864 + 1)
