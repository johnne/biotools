#!/usr/bin/env python

from Bio.SeqIO import parse
from Bio.SeqIO import write as seqio_write

from argparse import ArgumentParser
import sys

def initialize(record):
    records = []
    if len(str(record.seq))>=100: records.append(record)
    taxid=(record.id).split(".")[0]
    bioac=(record.id).split(".")[1]
    return [records,taxid,bioac]

def write_records(records,f,outdir):
    with open(outdir+"/"+f, "w") as output_handle:
        seqio_write(records, output_handle, "fasta")

def read_records_from_filehandle(fh,outdir):
    taxid = bioac = ""
    records = []
    for i,record in enumerate(parse(fh,"fasta")):
        if i==0: 
            [records,taxid,bioac] = initialize(record)
            continue
        thistaxid=(record.id).split(".")[0]
        thisbioac=(record.id).split(".")[1]
        if thistaxid!=taxid or thisbioac!=bioac:
            sys.stderr.write("Writing "+taxid+"."+bioac+" to "+outdir+"\n")
            write_records(records,taxid+"."+bioac+".fasta",outdir)
            records = []
            taxid = thistaxid
            bioac = thisbioac
        else:
            if len(str(record.seq))<100: continue
            #records+=[">"+record.id,str(record.seq)]
            records.append(record)
            taxid = thistaxid
            bioac = thisbioac
        

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Input fasta file from proGenomes (the representative.contigs.fasta.gz is a good choice)")
    parser.add_argument("-o", "--outdir", default=".",
            help="Write files to outdir. Defaults to current directory")
    args = parser.parse_args()

    with open(args.infile) as fh:
        read_records_from_filehandle(fh,args.outdir)

if __name__ == '__main__': 
    main()
