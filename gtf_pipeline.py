"""
This is a script for Ruffus pipeline

Function:
    - split file by chrom
    - count the number of transcript per chrom
    - read all count files, calculate the average and write to file
"""

#@
import sys
import gzip
from ruffus import *
from cgatcore import pipeline as P

@split('test.gtf.gz', 'chr*.gtf.gz')

def split_chrom(infile, outfiles):
    #P.run('sort k1,1 k4,4n %(infile)')
    with gzip.open(infile, 'rt') as inf:
        chrom = ''
        for line in inf:
            if line.split()[0] != chrom:
                chrom = line.split()[0]
                print(chrom)
                if chrom == '': outf.close()
                outf = gzip.open(chrom+'.gtf.gz', 'wt')
            outf.write(line)

#split_chrom('test.gtf.gz', '')

if __name__ == '__main__':
    sys.exit(P.main(sys.argv))
