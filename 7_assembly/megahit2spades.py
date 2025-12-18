#!/usr/bin/env python

import os
import sys
import argparse as ap
from Bio import SeqIO

def read_params(args):
    parser = ap.ArgumentParser(description='Convert from MEGAHIT to SPAdes format')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=sys.stdin, type=str, help="the input MEGAHIT file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str, help="the output SPAdes file [stdout if not present]")
    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    f = []
    fid = open(par['inp_f'], "r")
    for record in SeqIO.parse(fid,'fasta'):
        f.append(record)
    fid.close()

    if par['out_f']:
        fidout = open(par['out_f'],'w')
    else:
        fidout = sys.stdout

    l = [len(s) for s in f]
    t = sorted(range(len(l)), key=lambda s: l[s], reverse=True)
    c = 0
    for s in t:
        c = c+1
        f[s].id = 'NODE_' + str(c) + '_length_' + str(l[s]) + '_cov_' + f[s].description.split(' ')[2].split('=')[-1] + '_ID_' + f[s].id.split('_')[-1]
        f[s].description = ''
        SeqIO.write(f[s], fidout, "fasta")

    if par['out_f']:
        fidout.close()
