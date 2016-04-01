'''
This snakemake file prepares the input files and runs verifyBamID
to identify which cells are from which individual based on the RNA-seq
reads. Since in this experiment we know the identity of each cell, we
can determine how well this works.

The imputed Yoruba genotypes are in the directory:
/mnt/lustre/data/internal/genotypes/hg19/YRI. The genotypes are split
by chromosome and have the following naming scheme:
chr#.hg19.impute2.gz.  Each line of the file is a variant (includes
both SNPs and indels). The first column is for the chromosome, but
only contains --- as a placeholder. Thus the chromosome must be
obtained from the file name. The second column is the ID, the third is
the position, the fourth is the refence allele, and the fifth is the
alterate allele. There are three columns per sample. These contain the
probabilities that the sample is homozgyous reference, heterozygous,
or homozygous alternative. There are no column headers. The order of
the individuals is in the file YRI_samples.txt.

To submit:
nohup snakemake -kps verify-bam.py -j 96 --ri -c "qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log} -q blades.q" &
'''

import os

# Paths
GENOTYPES = "/mnt/lustre/data/internal/genotypes/hg19/YRI/"
BAM = "/mnt/gluster/home/jdblischak/ssd/bam-combined/"
EXONS = "/mnt/lustre/home/jdblischak/singleCellSeq/data/exons.saf"
DATA_DIR = 'data/'
LOG_DIR = 'log/'

for d in [DATA_DIR, LOG_DIR]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Settings
CHROM = [str(x) for x in range(1, 23)]

# Target rules
localrules: all, bam, vcf

# To-do: Compile all the output from verifyBamID into one file
rule all:
	input: "final.txt"

# Run verifyBamID on all the BAM files
rule bam:
	input: expand(DATA_DIR + "NA{IND}.r{REP}.{ROW}{COL}.bestSM", \
                  IND = ["19098", "19101", "19239"], \
                  REP = ["1", "2", "3"], \
                  ROW = ["A", "B", "C", "D", "E", "F", "G", "H"], \
                  COL = ["%02d"%(x) for x in range(1, 13)])

# Prepare the VCF files
rule vcf:
	input: expand(DATA_DIR + "chr{CHR}.hg19.exons.vcf", CHR = CHROM)

# Workflow ---------------------------------------------------------------------

# Convert IMPUTE2 genotypes to VCF format
# VCF format: http://samtools.github.io/hts-specs/VCFv4.3.pdf
rule convert_to_vcf:
	input: impute = GENOTYPES + 'chr{CHR}.hg19.impute2.gz',
           names = GENOTYPES + 'YRI_samples.txt'
	output: vcf = DATA_DIR + 'chr{CHR}.hg19.vcf'
	params: h_vmem = '2g', bigio = '0',
            name = lambda wildcards: 'convert_to_vcf.' + wildcards.CHR
	log: LOG_DIR
	run:
          import gzip
          import os

          out_handle = open(output.vcf, "w")

          # Get sample names
          names = []
          names_handle = open(input.names)
          for line in names_handle:
              cols = line.strip().split()
              names.append(cols[0])
          names_handle.close()

          header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
          out_handle.write(header)
          for n in names:
              out_handle.write("\t%s"%(n))
          out_handle.write("\n")

          f = input.impute
          chrom = os.path.basename(f).split(".")[0]
          impute_handle = gzip.open(f)
          for line in impute_handle:
              cols = line.decode().strip().split()
              snp_id = cols[1]
              position = cols[2]
              ref = cols[3]
              alt = cols[4]
              # Only include SNPs
              if len(ref) > 1 or len(alt) > 1:
                  continue
              # Output SNP information
              out_handle.write("%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\tGT"%(chrom, position, snp_id, ref, alt))
              # Determine genotypes
              geno_probs = [float(x) for x in cols[5:]]
              cutoff = 0.75
              for i in range(0, len(geno_probs), 3):
                  geno = "."
                  # Homozygous reference
                  if geno_probs[i] > cutoff:
                      geno = "0/0"
                  # Heterozygous
                  elif geno_probs[i + 1] > cutoff:
                      geno = "0/1"
                  # Homozygous reference
                  elif geno_probs[i + 2] > cutoff:
                      geno = "1/1"
                  out_handle.write("\t%s"%(geno))
              out_handle.write("\n")
          impute_handle.close()
          out_handle.close()


# Convert exons in SAF format to BED format. Duplicate exons are maintained.
rule convert_to_bed:
	input: saf = EXONS
	output: bed = DATA_DIR + "exons.bed"
	params: h_vmem = '2g', bigio = '0',
            name = "convert_to_bed"
	log: LOG_DIR
	run:
          saf = open(input.saf, "r")
          bed = open(output.bed, "w")
          # Discard header
          saf.readline()

          for line in saf:
              cols = line.strip().split("\t")
              id = cols[0]
              chr = cols[1]
              start = str(int(cols[2]) - 1)
              end = cols[3]
              strand = cols[4]
              entry = "%s\t%s\t%s\t%s\t%s\t%s\n"%(chr, start, end, id, 0, strand)
              bed.write(entry)

          saf.close()
          bed.close()

# Select only those SNPs in annotated exons of protein-coding genes
rule select_exonic_snps:
	input: vcf = DATA_DIR + "chr{CHR}.hg19.vcf",
           exons = DATA_DIR + "exons.bed"
	output: vcf = DATA_DIR + 'chr{CHR}.hg19.exons.vcf'
	params: h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'select_exonic_snps.' + wildcards.CHR
	log: LOG_DIR
	shell: "bedtools intersect -a {input.vcf} -b {input.exons} -u -header > {output.vcf}"

# Combine exonic SNPs into one file
rule combine_snps:
	input: vcf = expand(DATA_DIR + "chr{CHR}.hg19.exons.vcf", CHR = CHROM)
	output: vcf = DATA_DIR + 'snps.hg19.exons.vcf'
	params: h_vmem = '2g', bigio = '0',
            name = "combine_snps"
	log: LOG_DIR
	shell: "cat <(head -n 1 {input.vcf[0]}) <(cat {input.vcf} | grep -v CHROM) > {output.vcf}"

# Run verifyBamID to obtain the best individual match for the BAM file
rule verify_bam:
	input: vcf = DATA_DIR + "snps.hg19.exons.vcf",
           bam = BAM + "{IND}.{REP}.{ROW}{COL}.trim.sickle.sorted.combined.bam"
	output: DATA_DIR + "NA{IND}.r{REP}.{ROW}{COL}.bestSM"
	params: prefix = DATA_DIR + "NA{IND}.r{REP}.{ROW}{COL}",
            individual = "NA{IND}",
            h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'verify_bam.' + wildcards.IND
	log: LOG_DIR
	shell: "verifyBamID --vcf {input.vcf} --bam {input.bam} --best --ignoreRG --smID {params.individual} --out {params.prefix}"
