#!/usr/bin/env python2

# This script was modeled after the subsampling script I found on SeqAnswers
# for fasta/fastq files. I remodeled it for use with bam files by
# utilizing the pysam module.
# http://seqanswers.com/forums/showpost.php?p=60581&postcount=12

# usage:
#  python subsample-bam.py seed num outdir filename
#
#    seed - integer to seed random number generator
#    num - Number of reads to subsample
#    outdir - Path to directory to write file, is created if does not exist
#    filename - BAM file to be subsampled
#
# Should be used with Python 2
# Creates output file: outdir/basename.num.bam

################################################################################

import sys
if sys.version_info.major != 2:
    sys.stderr.write('failure\tMust use Python 2 to load pysam.\n')
    sys.exit(1)
import os
import random
import pysam

if __name__  == '__main__':
    
    seed = int(sys.argv[1])
    num = int(sys.argv[2])
    outdir = sys.argv[3]
    filename = sys.argv[4]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not os.path.exists(filename):
        sys.stderr.write('failure\t%s\tdoes not exist.\n'%(filename))
        sys.exit(1)
    
    bam = pysam.Samfile(filename, 'rb')

    total_count = bam.mapped + bam.unmapped
    if total_count < num:
        sys.stderr.write('failure\t%s\tcontains %d reads, which is less than requested: %d\n'%(
            filename, total_count, num))
        sys.exit(1)
    
    # Subsample num indices
    indices = list(xrange(total_count))
    random.seed(seed)
    random.shuffle(indices)
    subbed = indices[0:num]
    subbed.sort()
    
    output_name = outdir + '/' + os.path.basename(filename).rstrip('bam') + \
                  str(num) + '.bam'
    output = pysam.Samfile(output_name, 'wb', template = bam)
    index = 0
    list_index = 0
    for read in bam.fetch(until_eof = True):
        if index == subbed[list_index]:
            output.write(read)
            list_index += 1
            #sys.stderr.write(str(list_index) + '\n')
            if list_index == num:
                break
        index += 1
    output.close()
    
    pysam.index(output_name)
    # Test
    test = pysam.Samfile(output_name, 'rb')    
    if test.mapped + test.unmapped != num:
        sys.stderr.write('failure\t%s\nObserved: %d\n Expected: %d\n'%(
            output_name, test.mapped + test.unmapped, num))
        sys.exit(1)
    else:
        sys.stderr.write('success\t%s\n'%(test.filename))
    
    test.close()
