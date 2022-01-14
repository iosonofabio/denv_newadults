#!/usr/bin/env python
import os
import sys
import glob
import argparse
import subprocess as sp
import pandas as pd

datasets = [
    '20200313_4kids',
    '20200810_20kids',
    ]


data_fdn = '/home/fabio/projects/dengue_patients/sequencing_data/'
vdj_fdn = data_fdn+'vdj/'
summary_fdn = vdj_fdn+'summary/'

def find_samples(dataset):
    return os.listdir(f'{vdj_fdn}{dataset}')


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--dataset', choices=datasets, default=datasets[-1])
    runname = pa.parse_known_args()[0].dataset
    samplenames = find_samples(runname)

    os.makedirs(summary_fdn, exist_ok=True)

    tabs = []
    for sn in samplenames:
        print(sn)
        fdn_out = f'{vdj_fdn}/{runname}/{sn}/{sn}/outs/'
        if not os.path.isdir(fdn_out):
            print('VDJ output folder of cellranger not found, skip')
            continue

        fn_sum = f'{fdn_out}all_contig_annotations.csv'
        tab = pd.read_csv(fn_sum, sep=',')
        tab['patient_sample'] = sn
        tabs.append(tab)

    print('Concatenate outputs')
    tab_all = pd.concat(tabs)

    print('Write summary file')
    fn_sum = f'{summary_fdn}{runname}_vdj_all_contig_annotations.csv'
    tab_all.to_csv(fn_sum, index=False)
    print(f'Results writte in {fn_sum}')
