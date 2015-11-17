#!/usr/bin/env python
###################################
import sys,csv
from Bio.SeqIO import parse

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
            help="Input genbank file")
    parser.add_argument("-k", "--key", type=str, default="locus_tag",
            help="Unique gene key, defaults to 'locus_tag'")
    parser.add_argument("-f", "--feature", type=str, default="EC_number",
            help="Feature to extract and print for each unique key, defaults to 'EC_number'")
    parser.add_argument("-t", "--type", type=str, default="CDS",
            help="Type of feature to search, defaults to 'CDS'")
    args = parser.parse_args()

    if not args.infile: sys.exit(parser.print_help())

    from Bio.SeqIO import parse

    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')

    for record in parse(args.infile, 'genbank'):
        for feature in record.features:
            if not feature.type==args.type:continue
            qual = feature.qualifiers
            try: feats = qual[args.feature]
            except KeyError: continue
            try: key = qual[args.key]
            except KeyError: sys.exit("ERROR: Specified key "+args.key+" not found in input file\n")
            key = key[0]
            for feat in feats: houtcsv.writerow([key,feat])
    hout.close()

if __name__ == '__main__':
    main()
