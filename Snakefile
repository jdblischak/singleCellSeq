
"""

To run the full pipeline, submit the following line from within the
same directory as the Snakefile while on the head node (the paths to
the data files are relative to the Snakefile):

nohup snakemake -kp -j 96 --ri -c "qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log}" &

"""

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

localrules: test, qc

test_samples = ['lane1_Undetermined_L001_R1_001.fastq.gz',
                'lane1_Undetermined_L001_R1_007.fastq.gz',
                'lane2_Undetermined_L002_R1_031.fastq.gz']

rule test:
	input: [DATA_DIR + f.replace('fastq.gz', 'umi.bam') for f in test_samples]

rule qc:
	input: [DATA_DIR + f.replace('.fastq.gz', '_fastqc.zip') for f in test_samples]

###############################################################################
# Per-fastq processing: from raw reads to gene counts
###############################################################################

rule unzip:
	input: DATA_DIR + '{seq}.fastq.gz'
	output: DATA_DIR + '{seq}.fastq'
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
	output: DATA_DIR + '{seq}.trim.fastq'
	message: 'Trim UMIs from 5\' end of reads of sample {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'trim_umi.' + wildcards.seq
	log: LOG_DIR
	shell: '{UMITOOLS}umitools trim --end 5 {input} NNNNNGGG --verbose > {output}'

rule map:
	input: fastq = DATA_DIR + '{seq}.trim.fastq',
               genome = REF_GENOME + '.reads'
	output: DATA_DIR + '{seq}.bam'
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
