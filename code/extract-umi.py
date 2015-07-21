#!/usr/bin/env python

# Extract UMI from name of bam entries.
# Usage: extract-umi.py region outdir filename
#
# region - pull reads aligned to this region, e.g. ERCC-00130
# outdir - write files to this directory.
# filename - input bam file
#
################################################################################

import sys
if sys.version_info.major != 2:
    sys.stderr.write('failure\tMust use Python 2 to load pysam.\n')
    sys.exit(1)
import os
import pysam

if __name__  == '__main__':

    region = sys.argv[1]
    outdir = sys.argv[2]
    filename = sys.argv[3]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not os.path.exists(filename):
        sys.stderr.write('failure\t%s\tdoes not exist.\n'%(filename))
        sys.exit(1)
    
    bam = pysam.Samfile(filename, 'rb')

    output_name = outdir + '/' + os.path.basename(filename).rstrip('bam') + \
                  'umi.txt'
    output = open(output_name, 'w')
    output.write('chr\tstart\tumi\n')
    for read in bam.fetch(region, until_eof = True):
        chr = bam.getrname(read.tid)
        assert 'Read is aligned to specified region', chr == region
        start_pos = read.pos
        umi = read.qname.split(':')[-1].split('_')[-1]
        output.write(chr + '\t' + str(start_pos) + '\t' + umi + '\n')

    output.close()
    bam.close()
