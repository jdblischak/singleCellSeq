#!/usr/bin/env python

"""
Gathers all the counts contained in the summary files.

usage: gather-total-counts.py > total-counts.txt

Should be run from data directory.
"""

import glob
import sys

files = sorted(glob.glob("fastq/*count*")) + \
        sorted(glob.glob("trim/*count*")) + \
        sorted(glob.glob("sickle/*count*")) + \
        sorted(glob.glob("bam-processed/*count*")) + \
        sorted(glob.glob("bam-rmdup-umi/*count*")) + \
        sorted(glob.glob("counts/*summary"))

# files = glob.glob("fastq/*count*")[:2] + \
#         glob.glob("trim/*count*")[:2] + \
#         glob.glob("sickle/*count*")[:2] + \
#         glob.glob("bam-processed/*count*")[:2] + \
#         glob.glob("bam-rmdup-umi/*count*")[:2] + \
#         glob.glob("counts/*summary")[:2]

# print(len(files))

# Output header:
#   stage - the stage of the processing pipeline
#   sequences - reads or molecules
#   combined - TRUE if file is combination of all lanes for a given sample
sys.stdout.write("stage\tsequences\tcombined\tindividual\treplicate\twell\tindex\tlane\tflow_cell\tcounts\n")

# Function definitions ---------------------------------------------------------

def read_lines(f):
    """
    Input: Path to file
    Output: Lines of the file (list)
    """
    handle = open(f, "r")
    lines = handle.readlines()
    handle.close()
    return lines

def read_count_file(f):
    """
    Input: Path to file
    Output: The number contained in the file (str)
    Explanation: The count file only contains one number, the number of
                 sequences at that stage in the processing pipeline.
    """
    lines = read_lines(f)
    assert len(lines) == 1, "Count file has only one line"
    counts = lines[0].strip("\n")
    return counts

def read_featureCounts_summary(f):
    """
    Input: Path to file
    Output: The number of sequences labeled Assigned (str)
    Explanation: featureCounts outputs a file with the extension .summmary that
                 details the number of sequences per result category. The
                 category Assigned is for sequences that map unambiguously to a
                 gene.
    """
    assert f[-8:] == ".summary", \
           "featureCounts summary file has correct extension"
    lines = read_lines(f)
    assert len(lines) == 12, \
           "featureCounts summary file has 12 lines"
    assert lines[1].split("\t")[0] == "Assigned", \
           "The Assigned category is the first entry after the header"
    counts = lines[1].strip("\n").split("\t")[1]
    return counts

def check_combined(fname):
    """
    Input: Basename of file
    Output: TRUE if combined in fname, FALSE otherwise (str)
    Explanation: Files that include data from all lanes for a given sample
                 have combined in their name.
    """
    if "combined" in f:
        return "TRUE"
    else:
        return "FALSE"

# Process each file ------------------------------------------------------------

for f in files:
    # Set default values
    sequences = "reads"
    combined = "NA"
    index = "NA"
    lane = "NA"
    flow_cell = "NA"

    dir, fname = f.split("/")

    if dir == "fastq":
        stage = "raw"
        counts = read_count_file(f)
    elif dir == "trim":
        stage = "valid UMI"
        counts = read_count_file(f)
    elif dir == "sickle":
        stage = "quality trimmed"
        counts = read_count_file(f)
    elif dir == "bam-processed":
        stage = "mapped to genome"
        counts = read_count_file(f)
    elif dir == "bam-rmdup-umi":
        stage == "mapped to genome"
        sequences = "molecules"
        counts = read_count_file(f)
        combined = check_combined(fname)
    elif dir == "counts":
        stage = "mapped to exons"
        counts = read_featureCounts_summary(f)
        if "rmdup" in fname:
            sequences = "molecules"
        combined = check_combined(fname)

    # Get meta data from filename
    fname_parts = fname.split(".")
    if combined == "TRUE":
        individual, replicate, well = fname_parts[:3]
    else:
        individual, replicate, well, index, lane = fname_parts[:5]
        flow_cell = fname_parts[6]

    # Update variable names
    individual = "NA" + individual
    replicate = "r" + replicate

    sys.stdout.write(stage + "\t" + sequences + "\t" + combined + "\t" + \
                     individual + "\t" + replicate + "\t" + \
                     well + "\t" + index + "\t" + lane + "\t" + flow_cell + "\t" + \
                     counts + "\n")
