#!/usr/bin/env python
# John Blischak


'''
Sorts fastq files after demultiplexing with Casava.
'''

import sys
import os
import ipdb

################################################################################

def count_mismatch(x, y):
    num = 0
    assert len(x) == len(y), 'Strings must be same length.'
    assert isinstance(x, str) & isinstance(y, str), \
           'Input must be strings.'
    for i, j in zip(x, y):
        if i != j:
            num += 1
    return num

expected = open('expected_index.txt', 'r')

result = {}

for line in expected:
    exp = line.strip()
    result[exp] = [0] * 9

expected.close()


observed = open('index.fa', 'r')

for line in observed:
    if line[0] == '>':
        continue
    obs = line.strip()
    if obs == 'N' * 8:
        for index in result.keys():
            result[index][8] += 1
        continue
    for index in result.keys():
        mismatch = count_mismatch(index, obs)
        result[index][mismatch] += 1

observed.close()

out = open('mismatches.txt', 'w')
out.write('index\tmismatch\tcount\n')
for index in result.keys():
    for mismatch in range(9):
        out.write(index + '\t' + str(mismatch) + '\t' + \
                  str(result[index][mismatch]) + '\n')
out.close()
