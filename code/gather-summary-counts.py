#!/usr/bin/env python

# Gathers all the counts contained in the summary files.
# usage: gather-summary-counts.py > summary-counts.txt
# Should be run from data directory.

import glob
import sys

#files = glob.glob("counts/*combined*summary") # single cell molecules
files = glob.glob("counts/*.[A-H][0-1][0-9].*sorted.genecounts*summary") # single cell reads

# print(len(files))

sys.stdout.write("individual\treplicate\twell\trmdup\tAssigned\tUnassigned_Ambiguity\tUnassigned_MultiMapping\tUnassigned_NoFeatures\tUnassigned_Unmapped\tUnassigned_MappingQuality\tUnassigned_FragmentLength\tUnassigned_Chimera\tUnassigned_Secondary\tUnassigned_Nonjunction\tUnassigned_Duplicate\n")

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
    sys.stdout.write(individual + "\t" + replicate + "\t" + \
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
