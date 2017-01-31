#!/usr/bin/env python

import sys,pandas as pd
from argparse import ArgumentParser

def file_to_df(f):
    hin = open(f)
    r = []
    for line in hin:
        if line[0]=="#": continue
        line = line.rstrip()
        line = line.rsplit()
        items = [line[0],line[2],line[3],float(line[4]),float(line[5])]
        r.append(items)
    hin.close()
    df = pd.DataFrame(r)
    df.columns = ["query_seq","hmm_name","hmm_accession","evalue","score"]
    return df

def main():
    parser = ArgumentParser('''Simple parser that selects the highest scoring HMM per query. Users can set an e-value cutoff for filtering.''')
    parser.add_argument("-i","--infile",required=True,
            help="hmmer tabular output file")
    parser.add_argument("-e","--evalue",required=False,
            help="e-value cutoff to use")
    args = parser.parse_args()

    df = file_to_df(args.infile)

    if args.evalue: df = df.loc[df.evalue<args.evalue]

    df.sort_values("score",ascending=False,inplace=True)
    
    ann = df.groupby("query_seq").first()

    fh = sys.stdout
    ann.loc[:,["hmm_accession","evalue","score"]].to_csv(fh,sep="\t")

if __name__ == '__main__':
    main()
