#!/usr/bin/env python
import csv,sys
from argparse import ArgumentParser

def read_lengths(gff):
    hin = open(gff)
    hincsv = csv.reader(hin, delimiter = '\t')
    l = {}
    for row in hincsv:
        start = int(row[3])
        stop = int(row[4])
        length = stop-start+1
        g = row[-1].split(";")
        for item in g:
            if "gene_id" in item:
                gene = item.replace("gene_id ", "")
                break
        l[gene] = length
    hin.close()
    return l

def read_counts(files):
    c = {} ## Counts of genes in files
    g = {} ## Genes
    t = {} ## Total counts per file
    for f in files:
        c[f] = {}
        t[f] = 0
        hin = open(f)
        hincsv = csv.reader(hin, delimiter = '\t')
        for row in hincsv:
            c[f][row[0]] = int(row[1])
            t[f]+=int(row[1])
            g[row[0]] = ""
        hin.close()
    return [c,g,t]

def write_rpkms(c,g,t,l, raw=False):
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    houtcsv.writerow(["Gene"]+sorted(c.keys()))

    for gene, geneLength in l.iteritems():
        out = [gene]
        for f in sorted(c.keys()):
            totalNumReads = t[f]
            numReads = c[f][gene]
            if raw: rpkm = numReads
            else: 
                rpk = float(numReads) / (float(geneLength)/1000)
                rpkm = rpk / (float(totalNumReads)/1000000)
            out.append(rpkm)
        houtcsv.writerow(out)
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-c", "--countfile", type=str,
            help="Counts file (from e.g. htseq)")
    parser.add_argument("--counts", type=str,
            help="List of several count files to collate RPKMs for")
    parser.add_argument("-g", "--gff", type=str, required = True,
            help="GFF/GTF file used to define gene regions")
    parser.add_argument("-r", "--raw", action="store_true",
            help="Print raw read counts instead of computing RPKMs")

    args = parser.parse_args()

    if not args.countfile and not args.counts: 
        parser.print_help()
        sys.exit("\nERROR: Counts file(s) required\n")

    l = read_lengths(args.gff)

    if not args.counts: files = [args.countfile]
    else: files = [x.replace("\n","") for x in open(args.counts).readlines()]
    [c,g,t] = read_counts(files)

    write_rpkms(c,g,t,l, args.raw)

if __name__ == '__main__':
    main()

