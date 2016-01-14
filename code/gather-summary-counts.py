#!/usr/bin/env python

# Gathers all the counts contained in the summary files.
# usage: gather-summary-counts.py > summary-counts.txt
# Should be run from data directory.

import glob
import sys
import pandas as pd

# molecules + reads (for single cells only)
files = glob.glob("counts/*combined*summary") + \
        glob.glob("counts/*.[A-H][0-1][0-9].*sorted.genecounts*summary")

# print(len(files))

# First write to temporary file
tmp_file = open("/tmp/summary-counts-tmp.txt", "w")

tmp_file.write("individual\treplicate\twell\trmdup\tAssigned\tUnassigned_Ambiguity\tUnassigned_MultiMapping\tUnassigned_NoFeatures\tUnassigned_Unmapped\tUnassigned_MappingQuality\tUnassigned_FragmentLength\tUnassigned_Chimera\tUnassigned_Secondary\tUnassigned_Nonjunction\tUnassigned_Duplicate\n")

for f in files:
    dir, fname = f.split("/")
    d_counts = {}
    handle = open(f, "r")
    for line in handle:
        cols = line.strip("\n").split("\t")
        d_counts[cols[0]] = cols[1]
    handle.close()
    # Get meta data from filename
    fname_parts = fname.split(".")
    individual, replicate, well = fname_parts[:3]
    # Determine if sample is read or molecules counts, i.e. was
    # processed with umitools rmdup
    if "rmdup" in fname:
        rmdup = "molecules"
    else:
        rmdup = "reads"
    # Output meta data with the featureCounts summary data
    tmp_file.write(individual + "\t" + replicate + "\t" + \
                     well + "\t" + rmdup + "\t" + \
                     d_counts["Assigned"] + "\t" + \
                     d_counts["Unassigned_Ambiguity"] + "\t" + \
                     d_counts["Unassigned_MultiMapping"] + "\t" + \
                     d_counts["Unassigned_NoFeatures"] + "\t" + \
                     d_counts["Unassigned_Unmapped"] + "\t" + \
                     d_counts["Unassigned_MappingQuality"] + "\t" + \
                     d_counts["Unassigned_FragmentLength"] + "\t" + \
                     d_counts["Unassigned_Chimera"] + "\t" + \
                     d_counts["Unassigned_Secondary"] + "\t" + \
                     d_counts["Unassigned_Nonjunction"] + "\t" + \
                     d_counts["Unassigned_Duplicate"] + "\n")

tmp_file.close()

# Import tmp_file with pandas. Combine the reads per lane to reads per sample.
sum_counts_tmp = pd.read_table("/tmp/summary-counts-tmp.txt")

# Assert that there is the correct number of reads and molecules samples
num_reads_files = sum_counts_tmp.query('rmdup == "reads"').shape[0]
assert num_reads_files == 2592, \
    "Found %d reads files, expected 2592"%(num_reads_files)
num_molecules_files = sum_counts_tmp.query('rmdup == "molecules"').shape[0]
assert num_molecules_files == 864, \
    "Found %d molecules files, expected 864"%(num_molecules_files)

sum_counts = sum_counts_tmp.groupby(["individual", "replicate", "well", "rmdup"],
    as_index = False).sum()
sum_counts.to_csv(sys.stdout,
    sep = "\t", na_rep = "NA", index = False)
