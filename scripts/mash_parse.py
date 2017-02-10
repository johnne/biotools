#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd, sys

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Mash tabular file from running 'mash dist'")
    parser.add_argument("-e", "--evalue", default=1e-10, type=float,
            help="Threshold e-value. Filter hits with e-value above this threshold. Defaults to 1e-10")
    
    args = parser.parse_args()

    df = pd.read_csv(args.infile, header=None, sep="\t", names=["subject","query","distance","pvalue","shared"])
    
    df_filt = df.loc[df.pvalue<=args.evalue]
    df_filt.index = range(1,len(df_filt)+1)
    
    taxids = [x.split("-")[2] for x in df_filt.subject]
    taxdf = pd.DataFrame(taxids,columns=["taxid"])
    taxdf.index = df_filt.index
    
    df_filt = pd.merge(df_filt,taxdf,left_index=True,right_index=True)
    
    df_filt = df_filt.sort_values("distance")

    best_hits = df_filt.groupby("query").first()
    parsed_df = pd.DataFrame(columns=["taxid"])
    for q in df_filt["query"].unique(): 
        best = best_hits.loc[q]
        best_p = best["pvalue"]
        putative = df_filt.loc[(df_filt["query"]==q)&(df_filt["pvalue"]<=best_p*10),"taxid"]
        q_taxids = list(set([best["taxid"]]+list(putative.values)))
        q_taxdf = pd.DataFrame(index=[q]*len(q_taxids),columns=["taxid"],data=q_taxids)
        parsed_df = pd.concat([parsed_df,q_taxdf])

    parsed_df.to_csv(sys.stdout,sep="\t", header=False)

if __name__ == '__main__': main()
