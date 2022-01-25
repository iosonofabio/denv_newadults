'''
Author:	        Fabio Zanini
Date:	        2021/02/03
Description:	Call SNPs for the demux Male/Female

Needs conda environment bayes (freebayes + python>=3.7)
'''
import os
import sys
import argparse
import pandas as pd
import subprocess as sp


runsamples = {
    '20211218_4adults_duplexed': ['D_5194_5110', 'H_3008_3024'],
}

sampleconversion_10X = {
        'D_5194_5110': 'D_5194_5110_CKDL210026538-1a-SI_GA_C6_H2TW7DSX3',
        'H_3008_3024': 'H_3008_3024_CKDL210026537-1a-SI_GA_C5_H2TW7DSX3',
}


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--run', choices=tuple(runsamples.keys()), required=True)
    pa.add_argument('--samples', default=None)
    pa.add_argument('--steps', required=True)
    pa.add_argument('--dry', action='store_true')
    args = pa.parse_args()

    run = args.run

    if args.samples is not None:
        samples = args.samples.split(',')
    else:
        samples = runsamples[run]
    steps = args.steps.split(',')
    
    fdn_data = os.path.abspath('../../sequencing_data/')
    fdn_transcriptome = '/home/fabio/software/cellranger_reference_data/refdata-gex-GRCh38-2020-A'
 
    for sample in samples:
        print(sample)
        fdn_in = f'{fdn_data}/raw/{run}/raw_data/{sample}'
        fdn_out = f'{fdn_data}/bamfiles/{run}/{sample}'
        os.makedirs(fdn_out, exist_ok=True)

        if 'map' in steps:
            sample10X = sampleconversion_10X.get(sample, sample)
            cwd = os.path.abspath(os.curdir)
            call = 'cellranger-4.0.0 count '
            if args.dry:
                call += '--dry '
            try:
                call += (
                    '--nosecondary '
                    f'--id {sample} '
                    f'--transcriptome {fdn_transcriptome} '
                    f'--fastqs {fdn_in} '
                    f'--sample {sample10X} '
                    '--force-cells 8000 '
                    '--localcores 32 '
                    '--localmem 128 '
                    )
                print(call)
                os.chdir(fdn_out)
                #FIXME: just making sure
                #sp.run(call, shell=True, check=True)
            finally:
                os.chdir(cwd)

        if 'NtoD' in steps:
            print('Convert N CIGARs into D for freebayes - sigh!')
            import pysam

            fn_bam = f'{fdn_data}/bamfiles/{run}/{sample}/{sample}/outs/possorted_genome_bam.bam'
            fn_newbam = f'{fdn_data}/bamfiles/{run}/{sample}/possorted_genome_bam_NtoD.bam'
            os.makedirs(os.path.dirname(fn_newbam), exist_ok=True)
            with pysam.AlignmentFile(fn_bam, 'r') as f:
                # FIXME: just making sure
                if False:#not args.dry:
                    with pysam.AlignmentFile(fn_newbam, 'wb', template=f) as fout:
                        for ir, read in enumerate(f):
                            if ir % 100000 == 0:
                                print(ir, end='\r')
                            if read.cigarstring is None:
                                continue
                            if 'N' in read.cigarstring:
                                newcigar = []
                                for (code, length) in read.cigar:
                                    # 3 is N (intron), 2 is D (deletion)
                                    if code == 3:
                                        code = 2
                                    newcigar.append((code, length))
                                read.cigar = newcigar
                            fout.write(read)
                        print()

        if 'indexBAM' in steps:
            print('Make BAM index for freebayes')
            fn_newbam = f'{fdn_data}/bamfiles/{run}/{sample}/possorted_genome_bam_NtoD.bam'
            call = f'samtools index {fn_newbam}'
            print(call)
            if not args.dry:
                sp.run(call, shell=True, check=True)

        if 'freebayes' in steps:
            print('Call SNV via freebayes')
            vcf_fdn = f'{fdn_data}/vcf/{run}/{sample}'
            os.makedirs(vcf_fdn, exist_ok=True)
            fn_newbam = f'{fdn_data}/bamfiles/{run}/{sample}/possorted_genome_bam_NtoD.bam'
            fn_snv = f'{fdn_data}/vcf/{run}/{sample}/freebayes_snv.vcf'
            # TODO: adapt these files!
            fn_freebayes_genome = f'{fdn_transcriptome}/fasta/genome.fa'
            fn_freebayes_regions = f'{fdn_data}/vcf/genome_regions_freebayes.fa.fai'
            call = ('freebayes-parallel '
                    f'{fn_freebayes_regions} '
                    '32 ' # n CPUS
                    f'-f {fn_freebayes_genome} '
                    '-iXu ' # ignore indels
                    '-C 2 -q 1 ' # 2 reads, quality 1
                    f'{fn_newbam} '
                    f'> {fn_snv}'
                    )
            print(call)
            if not args.dry:
                sp.run(call, shell=True, check=True)

        if 'count_alleles' in steps:
            fn_bam = f'{fdn_data}/bamfiles/{run}/{sample}/{sample}/outs/possorted_genome_bam.bam'
            fn_barcodes = f'{fdn_data}/bamfiles/{run}/{sample}/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
            fn_snv = f'{fdn_data}/vcf/{run}/{sample}/freebayes_snv.vcf'
            fn_scSplit_bin = f'../../software/scSplit'
            fn_common_snps = f'../../software/scSplit_common_snvs_hg38.tar.gz'
            fdn_alleles = f'{fdn_data}/vcf/allele_counts_scSplit/{run}/{sample}'
            os.makedirs(fdn_alleles, exist_ok=True)
            call = (f'{fn_scSplit_bin} count '
                    f'-v {fn_snv} '
                    f'-i {fn_bam} '
                    f'-c {fn_common_snps} '
                    '-r ref_filtered.csv '
                    '-a alt_filtered.csv '
                    f'-o {fdn_alleles} '
                    f'-b {fn_barcodes}'
                    )
            print(call)
            if not args.dry:
                sp.run(call, shell=True, check=True)

            #import vcf
            #from collections import Counter
            #pos = Counter()
            #with open(fn_snv, 'r') as f:
            #    vcffile = vcf.Reader(f)
            #    for ir, record in enumerate(vcffile):
            #        pos[record.CHROM+':'+str(record.POS)] += record.samples[0].data.DP
            #        if ir == 10000:
            #            break

        if 'demux' in steps:
            fn_scSplit_bin = f'../../software/scSplit'
            fdn_alleles = f'{fdn_data}/vcf/allele_counts_scSplit/{run}/{sample}'
            fn_ref = f'{fdn_alleles}/ref_filtered.csv'
            fn_alt = f'{fdn_alleles}/alt_filtered.csv'
            fdn_demux = f'{fdn_data}/vcf/snv_demux/{run}/{sample}'
            os.makedirs(fdn_demux, exist_ok=True)
            call = (f'{fn_scSplit_bin} run '
                    f'-r {fn_ref} '
                    f'-a {fn_alt} '
                    '-n 2 '
                    f'-o {fdn_demux} '
                    )
            print(call)
            if not args.dry:
                sp.run(call, shell=True, check=True)


