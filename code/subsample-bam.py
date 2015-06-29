# John Blischak
# 27 Sep 2012

# MACS gives WARNINGS if the read counts between the treatment and control
# samples are drastically different because this could skew the FDR
# calculations (eg if the control has many more reads, it is more likely
# to have a peak by chance). Using the parameter --to-small to scale
# the larger data set towards the smaller one does not remove this
# warning. To overcome this, this script subsamples from the larger bam
# file the exact number of reads of the smaller sample.

# This script was modeled after the subsampling script I found on SeqAnswers
# for fasta/fastq files. I remodeled it for use with bam files by
# utilizing the pysam module.
# http://seqanswers.com/forums/showpost.php?p=60581&postcount=12

################################################################################

import pysam
import random
import sys
import os

if __name__  == '__main__':
    
    random.seed()
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    new_loc = os.getcwd()
    # For testing purposes:
    # file1 = r'/mnt/lustre/home/jdblischak/tag/sequences/JUNB/19238_C/seq.quality.sorted.nodup.bam'
    # file2 = r'/mnt/lustre/home/jdblischak/tag/sequences/JUNB/18358_C/seq.quality.sorted.nodup.bam'
    # For chimp:
    # file2 = r'/mnt/lustre/data/share/HCR_Chipseq/Mapped/Input/Input_C_Pooled36_111206.quality.sort.nodup.bam'
    # For human:
    # file2 = '/mnt/lustre/data/share/HCR_Chipseq/Mapped/Input/Input_H_Pooled36_111206.quality.sort.nodup.bam'
    # new_loc = r'/mnt/lustre/home/jdblischak/tag/sequences/JUNB/19238_C/'

    sys.stderr.write('Reading in bam files...\n')
    sam1 = pysam.Samfile(file1, 'rb')
    sam2 = pysam.Samfile(file2, 'rb')
    
    if sam1.mapped > sam2.mapped:
         larger = sam1
         smaller = sam2
    else:
         larger = sam2
         smaller = sam1
    sys.stderr.write('%s has more reads than %s.\n%d > %d\n'%(larger.filename,
                     smaller.filename, int(larger.mapped), int(smaller.mapped)))

    total = int(smaller.mapped)
    indices = list(xrange(larger.mapped))
    random.shuffle(indices)
    subbed = indices[0:total]
    subbed.sort()
    
    sys.stderr.write('Outputting subsampled bam file...\n')
    output_name = new_loc + '/' + larger.filename.split('/')[-1].strip('bam') + 'sub.' + str(smaller.mapped) + '.bam'
    output = pysam.Samfile(output_name, 'wb', template = larger)
    index = 0
    list_index = 0
    for read in larger.fetch():
        if index == subbed[list_index]:
            output.write(read)
            list_index += 1
            sys.stderr.write(str(list_index) + '\n')
            if list_index == total:
                break
        index += 1
    output.close()
    
    pysam.index(output_name)
    output = pysam.Samfile(output_name, 'rb')    
    if smaller.mapped != output.mapped:
        sys.stderr.write('Mistake. Reads in smaller dataset: %d. Reads in subsampled dataset: %d\n'%(smaller.mapped, output.mapped))
    else:
        sys.stderr.write('Finished successfully.\nWrote %s\n'%(output.filename))
    
    sam1.close()
    sam2.close()
    output.close()
    
         
    
    
