
"""

To run the full pipeline, submit the following line from within the
same directory as the Snakefile while on the head node (the paths to
the data files are relative to the Snakefile):

nohup snakemake -kp -j 96 --ri -c "qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log}" &

"""

from Bio import SeqIO

###############################################################################
# SOFTWARE
###############################################################################
FASTQC='/mnt/lustre/data/tools/FastQC/'
UMITOOLS='/mnt/lustre/home/jdblischak/programs/virtualenv-1.11/py2.7/bin/'
SAMTOOLS = '/usr/local/bin/'
SUBREAD = '/home/jdblischak/src/subread-1.4.4-Linux-x86_64/bin/'

###############################################################################
# Data
###############################################################################
DATA_DIR = 'data/seqs/'
LOG_DIR = 'log/'
REF_GENOME = 'data/genome/combined' # prefix only

###############################################################################
# Target rules
###############################################################################

localrules: all, test, qc

samples = glob_wildcards(DATA_DIR + '{seq}.fastq.gz')

test_samples = ['lane1_Undetermined_L001_R1_001.fastq.gz',
                'lane1_Undetermined_L001_R1_007.fastq.gz',
                'lane2_Undetermined_L002_R1_031.fastq.gz']

rule all:
	input: expand(DATA_DIR + '{seq}.counts.txt', seq = samples.seq)

rule test:
	input: [DATA_DIR + f.replace('fastq.gz', 'umi.bam') for f in test_samples]

rule qc:
	input: [DATA_DIR + f.replace('.fastq.gz', '_fastqc.zip') for f in test_samples]

###############################################################################
# Per-fastq processing: from raw reads to gene counts
###############################################################################

rule unzip:
	input: DATA_DIR + '{seq}.fastq.gz'
	output: temp(DATA_DIR + '{seq}.fastq')
	message: 'Unzipping sample {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'unzip.' + wildcards.seq
	log: LOG_DIR
	shell: 'zcat {input} > {output}'

rule fastqc:
	input: DATA_DIR + '{seq}.fastq'
	output: DATA_DIR + '{seq}_fastqc.zip'
	message: 'Running FastQC on sample {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'fastqc.' + wildcards.seq
	log: LOG_DIR
	shell: '{FASTQC}fastqc {input}'

rule trim_umi:
	input: DATA_DIR + '{seq}.fastq'
	output: temp(DATA_DIR + '{seq}.trim.fastq')
	message: 'Trim UMIs from 5\' end of reads of sample {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'trim_umi.' + wildcards.seq
	log: LOG_DIR
	shell: '{UMITOOLS}umitools trim --end 5 {input} NNNNNGGG --verbose > {output}'

rule map:
	input: fastq = DATA_DIR + '{seq}.trim.fastq',
               genome = REF_GENOME + '.reads'
	output: temp(DATA_DIR + '{seq}.bam')
	message: 'Map reads of sample {input.fastq}'
	params: h_vmem = '12g', bigio = '1',
	        name = lambda wildcards: 'map.' + wildcards.seq
	log: LOG_DIR
	shell: '{SUBREAD}subread-align -i {REF_GENOME} -r {input.fastq} --BAMoutput > {output}'

rule sort_bam:
	input: DATA_DIR + '{seq}.bam'
	output: DATA_DIR + '{seq}.sorted.bam'
	message: 'Sort bam file {input}'
	params: h_vmem = '8g', bigio = '1',
	        name = lambda wildcards: 'sort_bam.' + wildcards.seq,
                prefix = lambda wildcards: DATA_DIR + wildcards.seq + '.sorted'
	log: LOG_DIR
	shell: '{SAMTOOLS}samtools sort {input} {params.prefix}'

rule index_bam:
	input: DATA_DIR + '{seq}.sorted.bam'
	output: DATA_DIR + '{seq}.sorted.bam.bai'
	message: 'Index sorted bam file {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'index_bam.' + wildcards.seq
	log: LOG_DIR
	shell: '{SAMTOOLS}samtools index {input}'

rule rmdup_umi:
	input: bam = DATA_DIR + '{seq}.sorted.bam',
               index = DATA_DIR + '{seq}.sorted.bam.bai',
	output: bam = DATA_DIR + '{seq}.umi.bam',
                bed = DATA_DIR + '{seq}.umi.diffs.bed',
	message: 'Remove reads with duplicated UMIs for {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'rmdup_umi.' + wildcards.seq
	log: LOG_DIR
	shell: '{UMITOOLS}umitools rmdup {input.bam} {output.bam} > {output.bed}'

rule download_exons:
	output: 'data/genome/exons.saf'
	message: 'Create SAF file of human exons.'
	params: h_vmem = '8g', bigio = '0',
	        name = 'download_exons'
	log: LOG_DIR
	shell: 'Rscript download_exons.R > {output}'

rule ercc_tab:
	input: 'data/genome/ERCC92.fa'
	output: 'data/genome/ERCC92.txt'
	message: 'Create tab-sep text file of ERCC sequences.'
	params: h_vmem = '8g', bigio = '0',
	        name = 'ercc_tab'
	log: LOG_DIR
	run:
          out = open(output[0], 'w')
          for seq_record in SeqIO.parse("data/genome/ERCC92.fa", 'fasta'):
              out.write('%s\t%s\t1\t%d\t+\n'%(seq_record.id, seq_record.id,
                                              len(seq_record)))
          out.close()

rule combine_features:
	input: exons = 'data/genome/exons.saf',
               ercc = 'data/genome/ERCC92.txt'
	output: 'data/genome/exons_ERCC92.saf'
	message: 'Combine human exons and ERCC sequences.'
	params: h_vmem = '8g', bigio = '0',
	        name = 'combine_features'
	log: LOG_DIR
	shell: 'cat {input.exons} {input.ercc} > {output}'

rule featureCounts:
	input: reads = 'data/seqs/{seq}.sorted.bam',
               umi = 'data/seqs/{seq}.umi.bam',
               anno = 'data/genome/exons_ERCC92.saf'
	output: counts = 'data/seqs/{seq}.counts.txt',
                summary = 'data/seqs/{seq}.counts.txt.summary'
	message: 'Counts number of reads per feature for {input.umi}.'
	params: h_vmem = '8g', bigio = '1',
	        name = lambda wildcards: 'featureCounts.' + wildcards.seq
	log: LOG_DIR
	shell: '{SUBREAD}featureCounts -a {input.anno} -F SAF -o {output.counts} {input.reads} {input.umi}'

