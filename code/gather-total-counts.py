#!/usr/bin/env python

# Gathers all the counts contained in the summary files.
# usage: gather-total-counts.py > total-counts.txt
# Should be run from data directory.

import glob
import sys

files = glob.glob("fastq/*count*") + \
        glob.glob("trim/*count*") + \
        glob.glob("sickle/*count*") + \
        glob.glob("bam-processed/*count*") + \
        glob.glob("bam-rmdup-umi/*count*") + \
        glob.glob("counts/*summary")

# print(len(files))

sys.stdout.write("stage\tsickle\tindividual\tbatch\twell\tindex\tlane\tflow_cell\tcounts\n")

for f in files:
    # Ignore files that include data from all lanes for a given sample
    if "combined" in f:
        continue
    # Get counts from f
    handle = open(f, "r")
    if "summary" in f:
        lines = handle.readlines()
        counts = lines[1].strip("\n").split("\t")[1]
        dir, fname = f.split("/")
        if "rmdup" in fname:
          stage = dir + "-rmdup"
        else:
          stage = dir
    else:
      counts = handle.readline().strip("\n")
      dir, fname = f.split("/")
      stage = dir
    handle.close()
    # Get meta data from filename
    fname_parts = fname.split(".")
    individual, batch, well, index, lane = fname_parts[:5]
    flow_cell = fname_parts[6]
    # Determine if sample was quality trimmed with sickle
    if "sickle" in fname:
        sickle = "quality-trimmed"
    else:
        sickle = "not-quality-trimmed"
    # Samples in the sickle subdirectory are in the stage trim because they were just
    # processed with `umitools trim`.
    if stage == "sickle":
        stage = "trim"
    sys.stdout.write(stage + "\t" + sickle + "\t" + individual + "\t" + batch + "\t" + \
                     well + "\t" + index + "\t" + lane + "\t" + flow_cell + "\t" + \
                     counts + "\n")
