#!/usr/bin/env python

# Gathers all the counts per gene.
# usage: gather-gene-counts.py > gene-counts.txt
# Should be run from data directory.

import glob
import sys

files = glob.glob("counts/*genecounts.txt")

# print(len(files))

# Create header
sys.stdout.write("individual\tbatch\twell\tindex\tlane\tflow_cell\tsickle\trmdup")
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
    sys.stdout.write("\t" + gene)
handle.close()
sys.stdout.write("\n")

# Process each file
for f in files:
    # Get meta data from filename
    dir, fname = f.split("/")
    fname_parts = fname.split(".")
    individual, batch, well, index, lane = fname_parts[:5]
    flow_cell = fname_parts[6]
    # Determine if sample was quality trimmed with sickle
    if "sickle" in fname:
        sickle = "quality-trimmed"
    else:
        sickle = "not-quality-trimmed"
    # Determine if sample has had duplicate UMI-read pairs removed
    if "rmdup" in fname:
        rmdup = "molecules"
    else:
        rmdup = "reads"
    sys.stdout.write(individual + "\t" + batch + "\t" + \
                     well + "\t" + index + "\t" + lane + "\t" + flow_cell + "\t" + \
                     sickle + "\t" + rmdup)

    # Get counts from f
    g = 0 # iterator for indexing gene names
    handle = open(f, "r")
    for line in handle:
        if line[0] == "#" or line[:6] == "Geneid":
            continue
        cols = line.strip("\n").split("\t")
        gene = cols[0]
        assert gene == gene_list[g], "Gene names are not in correct order in file: %s"%(f)
        g += 1
        counts = cols[6]
        sys.stdout.write("\t" + counts)
    handle.close()
    sys.stdout.write("\n")
