'''
Snakemake pipeline for running subsample analyses.

To submit:
nohup snakemake -kps $ssc/code/subsample-pipeline.py -j 3600 --ri -c "qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log} -q blades.q" &
'''

import os

# Paths ------------------------------------------------------------------------
LOG_DIR = "log/"
COUNTS_MATRIX = "counts-matrix/"
FULL_DATA = "/mnt/lustre/home/jdblischak/singleCellSeq/data/"
GOOD_CELLS = "/mnt/lustre/home/jdblischak/singleCellSeq/data/quality-single-cells.txt"
KEEP_GENES = "/mnt/lustre/home/jdblischak/singleCellSeq/data/genes-pass-filter.txt"
OUTPUT_DIR = "output/"

for d in [LOG_DIR, OUTPUT_DIR]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Variables --------------------------------------------------------------------
individuals = ["NA19098", "NA19101", "NA19239"]
depths = [250000, 500000, 750000, 1000000, 1250000, 1500000]
cells = [5, 10, 15, 20, 25, 50, 75, 100, 125, 150]
seeds = range(1, 10 + 1)
types = ["reads", "molecules"]
quantiles = [0.5]

# Targets ----------------------------------------------------------------------
localrules: all, submit_subsampler, submit_subsampler_test

rule all:
	input: "subsampling-results.txt"

rule submit_subsampler:
    input: expand(OUTPUT_DIR + "{TYPE}-{IND}-{CELLS}-{SEED}-{DEPTH}.txt", \
           IND = individuals, \
           DEPTH = depths, \
           CELLS = cells, \
           SEED = seeds, \
           TYPE = types)

rule submit_subsampler_test:
    input: expand(OUTPUT_DIR + "{TYPE}-{IND}-{CELLS}-{SEED}-{DEPTH}.txt", \
           IND = individuals[0:2], \
           DEPTH = depths[0:2], \
           CELLS = cells[0:2], \
           SEED = seeds[0:2], \
           TYPE = types)

# Pipeline ---------------------------------------------------------------------
rule subsampler:
	input: single_sub = COUNTS_MATRIX + "{DEPTH}-{TYPE}-raw-single-per-sample.txt",
           single_full = FULL_DATA + "molecules-raw-single-per-sample.txt",
           ercc = FULL_DATA + "expected-ercc-molecules.txt",
           bulk = FULL_DATA + "reads-raw-bulk-per-sample.txt",
           good_cells = GOOD_CELLS,
           keep_genes = KEEP_GENES
	output: OUTPUT_DIR + "{TYPE}-{IND}-{CELLS}-{SEED}-{DEPTH}.txt"
	params: cells = "{CELLS}", seed = "{SEED}",
            ind = "{IND}",
            h_vmem = '2g', bigio = '0',
            name = lambda wildcards: 'sub-' + "-".join([wildcards.IND,
                                                       wildcards.CELLS,
                                                       wildcards.SEED,
                                                       wildcards.DEPTH])
	log: LOG_DIR
	shell: "\
subsampler.R {params.cells} {params.seed} {input.single_sub} \
  {input.single_full} {input.ercc} {input.bulk} \
  --individual={params.ind} --good_cells={input.good_cells} \
  --keep_genes={input.keep_genes} -d -o {output}"

rule gather_subsample_results:
	input: expand(OUTPUT_DIR + "{TYPE}-{IND}-{CELLS}-{SEED}-{DEPTH}.txt", \
           IND = individuals, \
           DEPTH = depths, \
           CELLS = cells, \
           SEED = seeds, \
           TYPE = types)
	output: "subsampling-results.txt"
	params: h_vmem = '2g', bigio = '0',
            name = "gather-subsampling-results"
	log: LOG_DIR
	run:
           import os

           # Grab the header from the first file
           first = open(input[0])
           header = first.readline()
           first.close()

           # Create output file
           out = open(output[0], "w")
           out.write("type\tdepth\t" + header)

           # Write the contents of each input file to output file
           for fname in input:
               fname_parts = os.path.basename(fname).rstrip(".txt").split("-")
               type = fname_parts[0]
               depth = fname_parts[-1]
               f = open(fname, "r")
               f_header = f.readline()
               assert "File %s has the correct header"%(fname), f_header == header
               for line in f:
                   out.write(type + "\t" + depth + "\t" + line)
               f.close()

           out.close()

