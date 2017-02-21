#!/usr/bin/env python

from Bio.SeqIO import parse
from Bio.SeqIO import write as seqio_write

from argparse import ArgumentParser
import sys

def write_seqs(records,f,outdir):
    with open(outdir+"/"+f, "w") as output_handle:
        seqio_write(records, output_handle, "fasta")

def read_records_from_filehandle(fh,outdir):
    taxid = bioac = ""
    records = []
    for record in parse(fh,"fasta"):
        newtaxid=(record.id).split(".")[0]
        newbioac=(record.id).split(".")[1]
        if newtaxid!=taxid or newbioac!=bioac and len(records)>0:
            sys.stderr.write("Writing "+taxid+"."+bioac+" to "+outdir+"\n")
            write_records(records,taxid+"."+bioac+".fasta",outdir)
            records = []
            taxid = newtaxid
            bioac = newbioac

        if len(str(record.seq))<100: continue
        #records+=[">"+record.id,str(record.seq)]
        records.append(record)
        

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Input fasta file from proGenomes (the representative.contigs.fasta.gz is a good choice)")
    parser.add_argument("-o", "--outdir", default=".",
            help="Write files to outdir. Defaults to current directory")
    args = parser.parse_args()

    with open(args.infile) as fh:
        read_records_from_filehandle(fh)
