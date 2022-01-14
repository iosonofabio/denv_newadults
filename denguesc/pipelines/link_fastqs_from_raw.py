import os
import sys
import glob
import argparse


datasets = [
    '20200810_20kids',
    ]


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--dataset', required=True, choices=datasets)
    args = pa.parse_args()

    dataset = args.dataset

    root_data_fdn = '/home/fabio/projects/dengue_patients/sequencing_data/'
    raw_fdn = f'{root_data_fdn}raw/{dataset}/'
    fastq_fdn = f'{root_data_fdn}fastq/{dataset}/'

    raw_fns = set(glob.glob(f'{raw_fdn}**/*.fastq.gz', recursive=True))

    fastq_fdns = os.listdir(fastq_fdn)
    print(fastq_fdns)

    missing = []
    for sn in fastq_fdns:
        fns = []
        for fn in raw_fns:
            if os.path.basename(fn).startswith(sn):
                fns.append(fn)
        for fn in fns:
            raw_fns.remove(fn)
        if len(fns) == 0:
            missing.append(sn)
            continue

        # Read 1 and 2 and lanes
        fn1, fn2 = [], []
        for fn in fns:
            if '_R1_' in fn:
                fn1.append(fn)
            elif '_R2_' in fn:
                fn2.append(fn)

        lanes = set()
        for fn in fn1:
            fn = os.path.basename(fn)
            i = fn.find('_L00')
            lane = fn[i+2:i+5]
            lanes.add(lane)

        if len(lanes) == 1:
            fdn_dest = f'{fastq_fdn}{sn}/'
            if len(fn1) != 1:
                print(fn1)
                raise RuntimeError('One lane but multiple read 1 fastq files')
            if len(fn2) != 1:
                print(fn2)
                raise RuntimeError('One lane but multiple read 2 fastq files')

            fnb1 = fdn_dest+os.path.basename(fn1[0])
            fnb2 = fdn_dest+os.path.basename(fn2[0])
            if not os.path.isfile(fnb1):
                os.link(fn1[0], fnb1)
            if not os.path.isfile(fnb2):
                os.link(fn2[0], fnb2)
        else:
            for lane in lanes:
                fdn_dest = f'{fastq_fdn}{sn}/L{lane}/'
                os.makedirs(fdn_dest, exist_ok=True)

                fn1l = [fn for fn in fn1 if f'_L{lane}_' in fn]
                fn2l = [fn for fn in fn2 if f'_L{lane}_' in fn]
                if len(fn1l) != 1:
                    print(fn1l)
                    raise RuntimeError('Multiple read 1 fastq files in the same lane')
                if len(fn2l) != 1:
                    print(fn2l)
                    raise RuntimeError('Multiple read 2 fastq files in the same lane')

                fnb1 = fdn_dest+os.path.basename(fn1[0])
                fnb2 = fdn_dest+os.path.basename(fn2[0])
                if not os.path.isfile(fnb1):
                    os.link(fn1[0], fnb1)
                if not os.path.isfile(fnb2):
                    os.link(fn2[0], fnb2)

            
