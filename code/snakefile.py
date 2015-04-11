

###############################################################################
FASTQ_DIR = 'fastq/'
BAM_DIR = 'bam/'
COUNTS_DIR = 'counts/'
LOG_DIR = 'log/'
REF_GENOME = 'genome/combined' # prefix only

###############################################################################
# Target rules
###############################################################################

localrules: all, qc

samples = glob_wildcards(FASTQ_DIR + '{seq}.fastq.gz')

rule all:
	input: expand(COUNTS_DIR + '{seq}.counts.txt', seq = samples.seq)

rule qc:
	input: expand(FASTQ_DIR + '{seq}_fastqc.zip', seq = samples.seq)

###############################################################################
# Per-fastq processing: from raw reads to gene counts
###############################################################################

rule unzip:
	input: FASTQ_DIR + '{seq}.fastq.gz'
	output: temp(FASTQ_DIR + '{seq}.fastq')
	message: 'Unzipping sample {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'unzip.' + wildcards.seq
	log: LOG_DIR
	shell: 'zcat {input} > {output}'

rule fastqc:
	input: FASTQ_DIR + '{seq}.fastq'
	output: FASTQ_DIR + '{seq}_fastqc.zip'
	message: 'Running FastQC on sample {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'fastqc.' + wildcards.seq
	log: LOG_DIR
	shell: 'fastqc {input}'

rule trim_umi:
	input: fastq = FASTQ_DIR + '{seq}.fastq',
               fastqc = FASTQ_DIR + '{seq}_fastqc.zip' # Not actually needed as input, but needs to be run
                                                       # as part of pipeline.
	output: fastq = temp(FASTQ_DIR + '{seq}.trim.fastq'),
                stats = FASTQ_DIR + '{seq}.filtered.umi.txt'
	message: 'Trim UMIs from 5\' end of reads of sample {input.fastq}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'trim_umi.' + wildcards.seq
	log: LOG_DIR
	shell: 'umitools trim --end 5 {input.fastq} NNNNNGGG --verbose --top 50 1> {output.fastq} 2> {output.stats}'

rule map:
	input: fastq = FASTQ_DIR + '{seq}.trim.fastq',
               genome = REF_GENOME + '.reads'
	output: temp(BAM_DIR + '{seq}.bam')
	message: 'Map reads of sample {input.fastq}'
	params: h_vmem = '12g', bigio = '1',
	        name = lambda wildcards: 'map.' + wildcards.seq
	log: LOG_DIR
	shell: 'subread-align -i {REF_GENOME} -r {input.fastq} --BAMoutput > {output}'

rule sort_bam:
	input: BAM_DIR + '{seq}.bam'
	output: BAM_DIR + '{seq}.sorted.bam'
	message: 'Sort bam file {input}'
	params: h_vmem = '16g', bigio = '1',
	        name = lambda wildcards: 'sort_bam.' + wildcards.seq,
                prefix = lambda wildcards: BAM_DIR + wildcards.seq + '.sorted'
	log: LOG_DIR
	shell: 'samtools sort {input} {params.prefix}'

rule index_bam:
	input: BAM_DIR + '{seq}.sorted.bam'
	output: BAM_DIR + '{seq}.sorted.bam.bai'
	message: 'Index sorted bam file {input}'
	params: h_vmem = '8g', bigio = '0',
	        name = lambda wildcards: 'index_bam.' + wildcards.seq
	log: LOG_DIR
	shell: 'samtools index {input}'

rule rmdup_umi:
	input: bam = BAM_DIR + '{seq}.sorted.bam',
               index = BAM_DIR + '{seq}.sorted.bam.bai',
	output: bam = BAM_DIR + '{seq}.umi.bam',
                bed = BAM_DIR + '{seq}.umi.diffs.bed',
	message: 'Remove reads with duplicated UMIs for {input.bam}'
	params: h_vmem = '16g', bigio = '0',
	        name =  lambda wildcards: 'rmdup_umi.' + wildcards.seq
	log: LOG_DIR
	shell: 'umitools rmdup {input.bam} {output.bam} > {output.bed}'

rule featureCounts:
	input: reads = BAM_DIR + '{seq}.sorted.bam',
               umi = BAM_DIR + '{seq}.umi.bam',
               anno = 'exons.saf'
	output: counts = COUNTS_DIR + '{seq}.counts.txt',
                summary = COUNTS_DIR + '{seq}.counts.txt.summary'
	message: 'Counts number of reads per feature for {input.umi}.'
	params: h_vmem = '8g', bigio = '1',
	        name = lambda wildcards: 'featureCounts.' + wildcards.seq
	log: LOG_DIR
	shell: 'featureCounts -a {input.anno} -F SAF -o {output.counts} {input.reads} {input.umi}'


