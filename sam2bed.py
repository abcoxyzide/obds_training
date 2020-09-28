# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:18:40 2020

@author: x5ei5

Write a Python script to convert the SAM file to a BED file

Problems:
    - POS is 1-based


What we want as output:
    - chr (col 3)
    - start pos (col 4 but 0 based)
    - end pos (col 3 + len(col10))
    - name (col 1)
    - score (col 5)
    - strand (.)

"""


import argparse

# # input file
# samfilepath = 'C:\\Users\\x5ei5\obds-wd\ERR1755082.test.sam' 
# # output file (doesn't exist before we created it)
# bedfilepath = 'C:\\Users\\x5ei5\obds-wd\ERR1755082.test.bed' 

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', dest='samfilepath', help='input file path')
parser.add_argument('--output', '-o', dest='bedfilepath', help='output file path')

args = parser.parse_args()

with open(args.samfilepath, 'r') as samfile: # open the input file (defined above)
    with open(args.bedfilepath, 'w') as bedfile: # create the output file
        for line in samfile: # read the input file line by line
            if line[0] == '@': # headers in SAM starts with @, so we skip it
                pass
            else:
                col = line.split() # parse the line into a list of fields / columns
                #     - chr (column 3)
                chrom = col[2]
                #     - start pos (column 4 but need to convert from 1-based SAM file to 0-based BED file)
                startpos = int(col[3]) - 1
                #     - end pos (column 4 + len(column 10))
                endpos = int(col[3]) + len(col[9])
                #     - name (column 1)
                name = col[0]
                #     - score (column 5)
                score = col[4]
                #     - strand (.)
                strand = '.'
                bedfile.write(f'{chrom}\t{startpos}\t{endpos}\t{name}\t{score}\t{strand}\n')






















