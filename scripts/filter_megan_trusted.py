#!/usr/bin/env python

import pandas as pd, sys
from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    parser.add_argument("-t","--trusted", required=True,
            help="List of trusted ORFs")
    parser.add_argument("-b", "--bed", required=True,
            help="BED file for assembly")
    parser.add_argument("-m", "--megan_file", required=True,
            helP="Megan taxid result file")
    args = parser.parse_args()
    
    trusted = pd.read_csv(args.trusted, header=None)
    bed = pd.read_csv(args.bed, header=None, sep="\t", index_col=0)
    megan = pd.read_csv(args.megan_file, header=None, index_col=0)

    contigs_with_trusted = list(set(bed[bed[3].isin(trusted[0])].index))

    ## Remove ORFs in the megan file that are: 1) not trusted and 2) on contigs with trusted ORFs
    notrust = set(bed[3]).difference(set(trusted[0]))
    orfs_contigs_with_trusted = bed.loc[contigs_with_trusted,3]
    to_remove = list(set(orfs_contigs_with_trusted).intersection(notrust))

    megan_filtered = megan.drop(to_remove, errors="ignore")
    megan_filtered.to_csv(sys.stdout)

if __name__ == '__main__':
    main()
