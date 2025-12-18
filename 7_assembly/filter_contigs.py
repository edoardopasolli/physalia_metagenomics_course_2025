#!/usr/bin/env python

import glob
import os
import sys
import argparse as ap
import pandas as pd
from Bio import SeqIO

def read_params(args):
    parser = ap.ArgumentParser(description='Filter contigs')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=sys.stdin, type=str, help="the input fasta file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str, help="the output file [stdout if not present]")
    arg( 'inp_f_filter', metavar='INPUT_FILE_FILTER', nargs='?', default=None, type=str, help="the input file with contigs to filter")
    arg( '-l','--minimum_read_length', default=1000, type=int, help="minimum read length")
    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    if par['inp_f_filter']:
        f_filter = pd.read_csv(par['inp_f_filter'], sep="\t")
        f_filter = f_filter["contig_name"].values

    f = []
    fid = open(par['inp_f'], "r")
    for record in SeqIO.parse(fid,'fasta'):
        if par['inp_f_filter']:
            if (len(record.seq) > par['minimum_read_length']) & (record.id not in f_filter):
                f.append(record)
        else:
            if (len(record.seq) > par['minimum_read_length']):
             f.append(record)
    fid.close()

    if par['out_f']:
        fid = open(par['out_f'],'w')
    else:
        fid = sys.stdout

    for s in range(len(f)):
        SeqIO.write(f[s], fid, "fasta")

    if par['out_f']:
        fid.close()
