#!/usr/bin/env python

import sys, gzip as gz, numpy as np, os
from argparse import ArgumentParser
from Bio.SeqIO import parse

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Fastq file")
    parser.add_argument("-m", "--maxseqs", required=False, type=int,
            help="Maximum number of sequences to use for average length calculation (optional)")
    args = parser.parse_args()

    l = []
    
    if ".gz" in args.infile: fh = gz.open(args.infile, 'rt')
    else: fh = open(args.infile, 'r')

    for i,record in enumerate(parse(fh, 'fastq'),start=1): 
        l.append(len(record.seq))
        if args.maxseqs and i>=args.maxseqs: break

    m = np.mean(l)
    name=os.path.basename(args.infile)
    print(name+'\t'+str(np.round(m,2)))

if __name__ == '__main__':
    main()
