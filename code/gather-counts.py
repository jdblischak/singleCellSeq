#!/usr/bin/env python

# Gathers all the counts contained in the summary files.
# usage: gather-counts.py > total-counts.txt
# Should be run from data directory.

import glob
import sys

files = glob.glob("trim/*count*") + \
        glob.glob("bam/*count*") + \
        glob.glob("bam-processed/*count*") + \
        glob.glob("bam-rmdup-umi/*count*") + \
        glob.glob("counts/*summary")

# print(len(files))

sys.stdout.write("stage\tindividual\tbatch\twell\tindex\tlane\tflow_cell\tcounts\n")

for f in files:
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
    sys.stdout.write(stage + "\t" + individual + "\t" + batch + "\t" + \
                     well + "\t" + index + "\t" + lane + "\t" + flow_cell + "\t" + \
                     counts + "\n")
