import os
import sys
import subprocess as sp


if __name__ == '__main__':

    print('Call change-o on 10X contigs for example patient')
    pname = '1-002'
    raw_fn = f'../../sequencing_data/vdj/20200810_20kids/{pname}/{pname}/outs/filtered_contig.fasta'
    input_fn = f'../../sequencing_data/vdj/changeo/input/patients/{pname}/filtered_contig.fasta'

    print('Copy file into new folder before igblast')
    sp.run(f'cp {raw_fn} {input_fn}', shell=True)

    print('Run igblast')
    igblast_fdn = '/home/fabio/projects/dengue/data/changeo/igblast'
    call = f'AssignGenes.py igblast -s {input_fn} -b {igblast_fdn} ' + \
            '--organism human --loci ig --format blast'
    print(call)
    sp.run(call, shell=True)

    print('Make db')
    call2 = f'MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s {input_fn} ' + \
             '-r IMGT_Human_*.fasta --10x filtered_contig_annotations.csv --extended'
    print(call2)

