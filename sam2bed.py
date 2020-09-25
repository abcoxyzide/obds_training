# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:18:40 2020

@author: x5ei5

Write a Python script to convert the SAM file to a BED file

Problems:
    - POS is 1-based in SAM file

What we want as output:
    - chr (col 3)
    - start pos (col 4 but 0 based)
    - end pos (col 3 + len(col10))
    - name (col 1)
    - score (col 5)
    - strand (.)

"""


import argparse, gzip, 
# import sys # for debugging purposes

# # input file
# samfilepath = 'C:\\Users\\x5ei5\obds-wd\ERR1755082.test.sam'
# # output file (doesn't exist before we created it)
# bedfilepath = 'C:\\Users\\x5ei5\obds-wd\ERR1755082.test.bed'

parser = argparse.ArgumentParser(description='Turn each read in SAM file into a BED format')
parser.add_argument('--input', '-i', dest='samfilepath', help='input file path')
parser.add_argument('--output', '-o', dest='bedfilepath', help='output file path')
parser.add_argument('--padding', '-p', dest='padding', type=int, default=0, help='amount of padding')
parser.add_argument('--minqual', '-q', dest='min_mapq', type=int, default=0, help='minimum MAPQ score')
parser.add_argument('--fragment', '-f', dest='fragment', action='store_true', help='a flag that outputs fragment instead of read coordinates')

args = parser.parse_args()
if args.bedfilepath[-3:] != '.gz': args.bedfilepath += '.gz' # append .gz extension if the supplied output has none

# sys.stdout.write(repr(args))

with open(args.samfilepath, 'r') as samfile: # open the input file (defined above)
    with gzip.open(args.bedfilepath, 'wt') as bedfile: # create the output file
        for line in samfile: # read the input file line by line
            if line[0] == '@': # headers in SAM starts with @, so we skip it
                pass
            else:
                col = line.split() # parse the line into a list of fields / columns
                # chr (column 3)
                chrom = col[2]
                if chrom == '*': # skip the rest of the loop (most importantly the bedfile.write) if chrom is unmapped
                    continue
                # score (column 5)
                score = int(col[4])
                if score > args.min_mapq: # only proceed if mapping score is good
                    # start pos (column 4 but need to convert from 1-based SAM file to 0-based BED file)
                    startpos = int(col[3]) - 1 - args.padding
                    # name (column 1)
                    name = col[0]
                    # strand (.)
                    strand = '.'
                    if args.fragment == True:
                        tlen = int(col[8])
                        if tlen > 0:
                            endpos = int(col[3]) + tlen
                        else:
                            continue
                    else:
                        endpos = int(col[3]) + len(col[9]) + args.padding

                    bedfile.write(f'{chrom}\t{startpos}\t{endpos}\t{name}\t{score}\t{strand}\n')
