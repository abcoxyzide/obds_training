"""
This is a script for Ruffus pipeline

Purpose:
- preprocessing of scRNAseq files

Flow:
1. cellranger

"""

import sys
import re
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P

params = P.get_parameters('pipeline_cellranger.yml')
samples = pd.read_csv("cellranger_samples.csv")
samples.set_index('name', inplace=True)
#print(samples)

@follows(mkdir("count"))
@transform('data/*/.sample', regex(r'data/(.+)/.sample'), r'count/\1/outs/filtered_feature_bc_matrix.h5')
def cellranger_count(infile, outfile):
    sampleid = re.search(r'data/(.+)/\.sample', infile).group(1)
    #print(sampleid)
    fastqs = samples['fastqs'][sampleid]
    cellnumber = samples['cells'][sampleid]
    chemistry = samples['chemistry'][sampleid]

    statement = """
    cellranger count
    --id=%(sampleid)s
    --transcriptome=%(cellrangercount_transcriptome)s
    --fastqs=%(fastqs)s
    --expect-cells=%(cellnumber)s
    --chemistry=%(chemistry)s
    --localcores=%(cellrangercount_threads)s
    --localmem=%(cellrangercount_memory)s
    > %(sampleid)s_standardout.log
    2> %(sampleid)s_standarderror.log
    && mv %(sampleid)s count/
    """
    totalmem=str(params["cellrangercount_memory"])+'G'
    P.run(statement, job_queue='all.q',
          job_threads=params["cellrangercount_threads"],
          job_total_memory=totalmem,
          job_condaenv='obds-py3')

if __name__ == '__main__':
    sys.exit(P.main(sys.argv))
