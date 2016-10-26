#!/usr/bin/env python
#####################

import sys, pandas as pd
from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    parser.add_argument("-c", "--cov", type=str, required=True,
            help="IMG coverage file with Average fold statistics per scaffold")
    parser.add_argument("-m", "--mapfile", type=str, required=True,
            help="Mapping file for scaffold names to JGI contigs")
    parser.add_argument("-g", "--gff", type=str, required=True,
            help="GFF file for contigs")

    args = parser.parse_args()

    nmap = pd.read_csv(args.mapfile, header=None, sep="\t", index_col=0)
    rename = nmap.to_dict()

    gff = pd.read_csv(args.gff, header=None, sep="\t", usecols=[0,8], names=["Contig","Gene"], index_col=0)
    tags = [x.split("=")[-1].rstrip(";") for x in gff.Gene]
    gff["Gene"] = tags

    ## Make cov table and contigs in the GFF match up
    cov = pd.read_csv(args.cov, header=0, sep="\t", index_col=0)
    cov.rename(index = rename[1], inplace=True)

    ## Transfer coverage info to locus tags
    tmp = pd.merge(gff,cov,left_index=True, right_index=True)
    gene_cov = tmp.loc[:,["Gene","Avg_fold"]]

    gene_cov.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == '__main__':
    main()
