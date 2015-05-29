#!/usr/bin/env python

# Gathers all the counts contained in the summary files.
# usage: gather-summary-counts.py > summary-counts.txt
# Should be run from data directory.

import glob
import sys

files = glob.glob("counts/*summary")

# print(len(files))

sys.stdout.write("individual\tbatch\twell\tindex\tlane\tflow_cell\trmdup\tsickle\tAssigned\tUnassigned_Ambiguity\tUnassigned_MultiMapping\tUnassigned_NoFeatures\tUnassigned_Unmapped\tUnassigned_MappingQuality\tUnassigned_FragementLength\tUnassigned_Chimera\tUnassigned_Secondary\n")

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
    individual, batch, well, index, lane = fname_parts[:5]
    flow_cell = fname_parts[6]
    # Determine if sample is read or molecules counts, i.e. was
    # processed with umitools rmdup
    if "rmdup" in fname:
        rmdup = "molecules"
    else:
        rmdup = "reads"
    # Determine if sample was quality trimmed with sickle
    if "sickle" in fname:
        sickle = "quality-trimmed"
    else:
        sickle = "not-quality-trimmed"
    # Output meta data with the featureCounts summary data
    sys.stdout.write(individual + "\t" + batch + "\t" + \
                     well + "\t" + index + "\t" + lane + "\t" + flow_cell + "\t" + \
                     rmdup + "\t" + sickle + "\t" + \
                     d_counts["Assigned"] + "\t" + \
                     d_counts["Unassigned_Ambiguity"] + "\t" + \
                     d_counts["Unassigned_MultiMapping"] + "\t" + \
                     d_counts["Unassigned_NoFeatures"] + "\t" + \
                     d_counts["Unassigned_Unmapped"] + "\t" + \
                     d_counts["Unassigned_MappingQuality"] + "\t" + \
                     d_counts["Unassigned_FragementLength"] + "\t" + \
                     d_counts["Unassigned_Chimera"] + "\t" + \
                     d_counts["Unassigned_Secondary"] + "\n")
