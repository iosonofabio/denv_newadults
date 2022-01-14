import os
import sys
import numpy as np
import pandas as pd
import pysam
from collections import Counter, defaultdict
import argparse


runs = [
    'ZY_10010001',
]


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--run', choices=runs, required=True)
    pa.add_argument('--nbarcodes', default=-1, type=int)
    pa.add_argument('--nreads', default=-1, type=int)
    pa.add_argument('--test', action='store_true')
    pa.add_argument('--buffersize', default=10000, type=int)
    args = pa.parse_args()

    fdn = f'../../sequencing_data/bamfiles/{run}/'
    if args.test:
        fnmuxed = fdn+'test/head_test.bam'
        fdndemux = 'test/demux/'
    else:
        fnmuxed = fdn+'possorted_genome_bam.bam'
        fdndemux = 'demux/'

    barcodes_fn = f'../../sequencing_data/fastq/{run}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    barcodes = pd.read_csv(barcodes_fn, sep='\t', compression='gzip', squeeze=True).values
    if args.nbarcodes != -1:
        barcodes = barcodes[:args.nbarcodes]

    def remove_all_demuxed_bamfiles():
        for fn in os.listdir(fdn+fdndemux):
            os.remove(fdn+fdndemux+fn)

    def write_demux_bam(cb, reads, nfiles, template):
        if nfiles == 0:
            fn = fdn+fdndemux+'{:}.bam'.format(cb)
        else:
            fn = fdn+fdndemux+'{:}_{:}.bam'.format(cb, nfiles+1)
        with pysam.AlignmentFile(fn, 'wb', template=template) as bamfile:
            for read in reads:
                bamfile.write(read)

    def merge_bamfiles(cb, nfiles):
        import subprocess as sp

        fns = [fdn+fdndemux+'{:}_{:}.bam'.format(cb, x+1) for x in range(nfiles)]
        fn_out = fdn+fdndemux+'{:}.bam'.format(cb)
        # The first file needs renaming
        os.rename(fn_out, fns[0])

        # Merge
        sp.run('samtools merge {:} '.format(fn_out)+' '.join(fns), shell=True)

        # Delete tmp files
        for fn in fns:
            os.remove(fn)

    def flush_buffer(buf, bufnfiles, bamfile):
        ncb = len(buf)
        for i, (cb, reads) in enumerate(buf.items()):
            print('Flushing {:} / {:}'.format(i+1, ncb), end='\r')
            write_demux_bam(cb, reads, bufnfiles[cb], bamfile)
            if bufnfiles[cb] > 0:
                merge_bamfiles(cb, bufnfiles[cb]+1)
        print()

    if args.test:
        print('Count reads')
        with pysam.AlignmentFile(fnmuxed, 'rb') as bamfile:
            nreads = sum(1 for read in bamfile)
        print('Total number of reads: {:}'.format(nreads))
    else:
        if args.nreads:
            nreads = args.nreads
        else:
            nreads = 2e9

    def demux_reads(fnmuxed, nreads):
        remove_all_demuxed_bamfiles()

        buf = defaultdict(list)
        bufn = Counter()
        bufnfiles = Counter()
        with pysam.AlignmentFile(fnmuxed, 'rb') as bamfile:
            for ir, read in enumerate(bamfile):
                if ((ir + 1) % 100000) == 0:
                    print('read {:} / {:}'.format(ir+1, nreads), end='\r')
                if ir == args.nreads:
                    break

                if not read.has_tag('CB'):
                    continue

                cb = read.get_tag('CB')
                if cb not in barcodes:
                    continue

                buf[cb].append(read)
                bufn[cb] += 1
                if bufn[cb] == args.buffersize:
                    write_demux_bam(cb, buf[cb], bufnfiles[cb], bamfile)
                    buf[cb] = []
                    bufn[cb] = 0
                    bufnfiles[cb] += 1
            print()

            print('Flushing buffer and merging output files')
            flush_buffer(buf, bufnfiles, bamfile)

    demux_reads(fnmuxed, nreads)

    print('Check what fraction of barcodes actually get a BAM file')
    bamfile_found = pd.Series(np.zeros(len(barcodes), bool), index=barcodes)
    for fn in os.listdir(fdn+fdndemux):
        cb = fn.split('.')[0]
        bamfile_found[cb] = True
    print('Fraction of barcode bamfiles: {:.2f}'.format(bamfile_found.mean()))
