"""
This is a script for Ruffus pipeline

Function:
    - split file by chrom
    - count the number of transcript per chrom
    - read all count files, calculate the average and write to file
"""

import sys
import gzip
from ruffus import *
from cgatcore import pipeline as P

@split('genes.gtf.gz', 'chr*.gtf.gz')
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

@transform('chr*.gtf.gz', suffix('.gtf.gz'), '.count')
def count_genes(infile, outfile):
    statement = "wc -l %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@merge(count_genes, 'all.average')
def average(infiles, outfile):
    counts = []
    for infile in infiles:
        with open(infile, 'r') as inf:
            counts.append(int(inf.read().split()[0]))
    average = sum(counts) / len(counts)
    with open(outfile, 'w') as outf:
        outf.write(f'The average transcript count is {average}\n')


if __name__ == '__main__':
    sys.exit(P.main(sys.argv))
