#!/usr/bin/env python

import sys, csv
from argparse import ArgumentParser

def Parse(infile, i,e,l):
    hin = open(infile, 'r')
    hits = {}
    for row in hin:
        row = row.rstrip()
        row = row.rsplit()
        [qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore] = row[0:12]
        if float(pident)<i or float(evalue)>=e or int(length)<l: continue
        sstart = int(sstart)-1
        send = int(send)-1
        if sstart>send: send,sstart = sstart,send
        hit = [sseqid,str(sstart),str(send),qseqid]
        try: hits[qseqid].append(hit)
        except KeyError: hits[qseqid] = [hit]
    hin.close()
    return hits

def Print(hits):
    for query,l in hits.iteritems():
        for item in l: print "\t".join(item)

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
            help="Blast tabular output")
    parser.add_argument("-e", "--evalue", type=float, default=1e-5,
            help="Maximum e-value. Defaults to 1e-5")
    parser.add_argument("-I", "--id", type=float, default=100.0,
            help="Minimum percent identity. Defaults to 100")
    parser.add_argument("-a", "--alignlen", type=int, default=60,
            help="Minimum alignment length. Defaults to 60")

    args = parser.parse_args()
    hits = Parse(args.infile, args.id,args.evalue,args.alignlen)
    Print(hits)

if __name__ == '__main__':
    main()
