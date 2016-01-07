#!/usr/bin/env python

# Gathers all the counts per gene.
# usage: gather-gene-counts.py [files] > gene-counts.txt
#   e.g. gather-gene-counts.py counts/*genecounts.txt > gene-counts.txt
# Should be run from data directory.

# Creates 6 files:
#
# 1. reads-raw-bulk-per-lane.txt
# 2. reads-raw-bulk-per-sample.txt
# 3. reads-raw-single-per-lane.txt
# 4. reads-raw-single-per-sample.txt
# 5. molecules-raw-single-per-lane.txt
# 6. molecules-raw-single-per-sample.txt

import glob
import sys
import pandas as pd

if len(sys.argv) > 1:
    files = sys.argv[1:]
else:
    files = glob.glob("counts/*genecounts.txt")

# print(len(files))

# Files 1, 3, 5, and 6 can be created directly from the data.

reads_raw_bulk_per_lane = open("reads-raw-bulk-per-lane.txt", "w")
reads_raw_single_per_lane = open("reads-raw-single-per-lane.txt", "w")
molecules_raw_single_per_lane = open("molecules-raw-single-per-lane.txt", "w")
molecules_raw_single_per_sample = open("molecules-raw-single-per-sample.txt", "w")

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
