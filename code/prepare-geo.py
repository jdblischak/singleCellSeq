#!/usr/bin/env python


import glob
import string
import re
import hashlib

fastq_all_full = glob.glob("fastq/*fastq.gz")
fastq_all = [f.split("/")[-1] for f in fastq_all_full]
fastq_all.sort()
  
individual = ["19098", "19101", "19239"]
replicate = ["1", "2", "3"]
well = []
for row in string.ascii_uppercase[:8]:
    for col in range(1, 13):
        well = well + ["%s%02d"%(row, col)]
well = well + ["bulk"]

################################################################################
# Samples
################################################################################

samples = open("geo-samples.txt", "w")

# Static columns
source_name = "LCL-derived iPSC"
organism = "Homo sapiens"
molecule = "polyA RNA"
description = ""

# Sample counter
i = 1

for ind in individual:
    for rep in replicate:
        for w in well:
            # meta data
            title = "-".join(["NA" + ind, "r" + rep, w])
            samples.write("\t".join([
                "Sample %d"%(i),
                title,
                source_name,
                organism,
                "NA" + ind,
                "r" + rep,
                w,
                molecule,
                description]))
            # Processed data files
            if w == "bulk":
                samples.write("\t" + "\t".join([
                "reads-raw-bulk-per-lane.txt",
                "reads-raw-bulk-per-sample.txt ",
                "",
                ""]))
            else:
                samples.write("\t" + "\t".join([
                "reads-raw-single-per-lane.txt",
                "reads-raw-single-per-sample.txt",
                "molecules-raw-single-per-lane.txt",
                "molecules-raw-single-per-sample.txt"]))                
            # Raw data files
            search_term = ".".join([ind, rep, w])
            fastq = [fname for fname in fastq_all \
                     if re.search(search_term, fname)]

            if w == "bulk":
                num_raw_bulk = len(fastq)
                assert num_raw_bulk == 8, \
                    "Need 8 raw bulk files, only %d found"%(num_raw_bulk)
                samples.write("\t" + "\t".join(fastq))
            else:
                num_raw_single = len(fastq)
                assert num_raw_single == 3, \
                    "Need 3 raw single files, only %d found"%(num_raw_single)
                samples.write("\t" + "\t".join(fastq + [""] * 5))

 
            # Update sample counter            
            i +=1

            samples.write("\n")
        

samples.close()

################################################################################
# Processed data files
################################################################################

processed = open("geo-processed-data-files.txt", "w")

# static column
file_type = "tab-delimited text"

processed_files = glob.glob("counts-matrix/*txt")
processed_files.sort()

for p in processed_files:
    p_md5 = hashlib.md5(open(p, 'r').read()).hexdigest()
    fname = p.split("/")[1]
    processed.write("\t".join([fname, file_type, p_md5]) + "\n")

processed.close()

################################################################################
# Raw files
################################################################################

raw = open("geo-raw-files.txt", "w")

# static columns
file_type = "fastq"
instrument_model = "HiSeq 2500"
read_length = "100"
single_or_paired_end = "single"

# For each fastq file, read the previously computed md5sum hash
# http://jdblischak.github.io/singleCellSeq/analysis/verify-md5sum.html#calculate-md5-checksums
for f in fastq_all:
    md5file = open("md5sum/" + \
                   f.replace("fastq.gz", "md5.txt"))
    md5sum = md5file.read().strip().split()[0]
    md5file.close()
    raw.write("\t".join([
        f,
        file_type,
        md5sum,
        instrument_model,
        read_length,
        single_or_paired_end]) + "\n")

raw.close()
