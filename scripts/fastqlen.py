#!/usr/bin/env python

import sys, gzip as gz, numpy as np
from argparse import ArgumentParser
from Bio.SeqIO import parse

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Fastq file")
    args = parser.parse_args()

    l = []
    
    if ".gz" in args.infile: fh = gz.open(infile, 'rt')
    else: fh = open(infile, 'r')

    for record in parse(fh, 'fastq'): l.append(len(record.seq))

    m = np.mean(l)

    print(args.infile+'\t'+str(np.round(m,2)))

if __name__ == '__main__':
    main()
