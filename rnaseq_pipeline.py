"""
This is a script for Ruffus pipeline

Purpose:
- process fastq RNAseq files into counts

Flow:
1. fastQC
2. multiQC
3.

"""

import sys
import gzip
from ruffus import *
from cgatcore import pipeline as P

params = P.get_parameters('rnaseq_pipeline.yml')

@follows(mkdir("fastqc"))
@transform('*fastq.gz', regex(r'(.*)\.fastq\.gz'), r'fastqc/\1_fastqc.html')
def fastqc(infile, outfile):
    statement = "fastqc --nogroup -o fastqc %(infile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@merge(fastqc, 'fastqc/multiqc_report.html')
def multiqc(infiles, outfile):
    statement = "multiqc -f -n %(outfile)s fastqc"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@follows(mkdir("bam"))
# @collate( input, filter, output, [extras,...] )
@collate('*.fastq.gz', regex(r'(.*)_[12](.*)\.fastq\.gz'), r'bam/\1\2.bam')
def mapping(infiles, outfile):
    read1, read2 = infiles
    statement = """
    hisat2 -x %(hisat2_index)s
    -1 %(read1)s
    -2 %(read2)s
    --rna-strandness %(hisat2_strandness)s
    --threads %(hisat2_threads)s
    --summary-file %(outfile)s.log |
    samtools sort - -o %(outfile)s
    """
    P.run(statement, job_queue='all.q', job_threads=params["hisat2_threads"], job_memory='2G', job_condaenv='obds-py3')

@transform(mapping, regex(r'(.*)\.bam'), r'\1.idxstats')
def idxstats(infile, outfile):
    statement = "samtools idxstats %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@transform(mapping, regex(r'(.*)\.bam'), r'\1.flagstat')
def flagstat(infile, outfile):
    statement = "samtools flagstat %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@merge(mapping, 'readcount.tsv')
def featcount(infiles, outfile):
    inf = ' '.join(infiles)
    statement = """
    featureCounts %(featureCounts_opts)s
    -a %(featureCounts_annot)s
    -T %(featureCounts_threads)s
    -s %(featureCounts_strandness)s
    %(inf)s
    -o %(outfile)s
    """
    P.run(statement, job_queue='all.q', job_threads=params["featureCounts_threads"], job_memory='2G', job_condaenv='obds-py3')

@merge([idxstats, flagstat, featcount], 'bam/multiqc_report_bam.html')
def bamultiqc(infiles, outfile):
    statement = "multiqc -f -n %(outfile)s bam"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')


if __name__ == '__main__':
    sys.exit(P.main(sys.argv))
