#!/usr/bin/env python

import pandas as pd, sys
from argparse import ArgumentParser
from skbio.tree import TreeNode, MissingNodeError

def load_ranks(nodes_file):
    ranks = {}
    with open(nodes_file) as nodes_f:
        for line in nodes_f:
            line = line.strip()
            data_fields = line.split('\t')
            ranks[int(data_fields[0])] = {'rank': data_fields[4]}
    return pd.DataFrame.from_dict(ranks, orient='index')

def assign_ancestors(taxids,ncbi_tree,ranks_df,val_cols):
    default_val_col_values = [None] * len(val_cols)
    unknown_val_col_values = ["Unknown"] * len(val_cols)
    taxdf = {}
    for taxid in taxids:
        ## Initialize row and add default values
        data_row = {'tax_id': taxid}
        data_row.update(zip(val_cols, default_val_col_values))
        ## Find the node in the ncbi tree structure
        try: node = ncbi_tree.find(str(taxid))
        ## If the node is missing, set "Unknown" values for all ranks and continue
        except MissingNodeError: 
            data_row.update(zip(val_cols,unknown_val_col_values))
            continue
        ## Find and add each ancestor, if the ancestor rank is to be stored
        for ancestor in node.ancestors():
            try: rank = ranks_df['rank'].ix[int(ancestor.name)]
            except KeyError: continue
            except TypeError: continue
            if rank in val_cols: data_row[rank] = ancestor.name
        ## Find and add the taxid itself to the data row
        try: taxid_rank = ranks_df['rank'].ix[int(taxid)]
        except KeyError: taxid_rank = "Unknown"
        if taxid_rank in val_cols: data_row[taxid_rank] = taxid      
        taxdf[taxid] = data_row
    return pd.DataFrame.from_dict(taxdf,orient='index')

def lca(result_df,taxdf,ncbi_map,val_cols):
    lcadf = {}
    unknown_val_col_values = ["Unknown"] * len(val_cols)    
    for contig_id in result_df.index.unique():
        ## Initialize the row for the contig
        lca_res = {}        
        ## Get all taxids assigned to this contig
        taxids = result_df.loc[contig_id]
        ## Handle case where there's only one taxid assigned
        if len(taxids)==1: taxids = [taxids['taxid']]
        else: taxids = taxids['taxid'].values
        ## Check to see that the taxids have ancestor assignments
        try: taxa = taxdf.loc[taxids]
        ## If not, add "Unknown" default values
        except KeyError: 
            lca_res.update(zip(val_cols,unknown_val_col_values))
            continue
        ## Set LCA index to the lowest rank in val_cols
        lca_index = len(val_cols)-1
        ## Iterate ranks (assumes they are sorted from highest to lowest)
        for i,rank in enumerate(val_cols):
            ## For each rank, get assigned taxids
            rank_ids = list(set(taxa[rank].values))
            ## If there's more than one taxid, the previous rank is the LCA
            if len(rank_ids)>1:
                lca_index = i-1
                break
            ## Otherwise, add the taxa name for this rank
            else: 
                try: lca_res[rank] = ncbi_map['name'].ix[int(rank_ids[0])]
                except KeyError: lca_res[rank] = rank_ids[0]
                except TypeError: lca_res[rank] = "" ## Handle case with missing values for this specific rank (*)
        if lca_index == -1: lca_name = ""
        else: lca_name = lca_res[val_cols[lca_index]]
        ## Assign names to ranks below the LCA
        if lca_index<len(val_cols)-1:
            for rank in val_cols[lca_index+1:]: lca_res[rank] = ("Unclassified."+lca_name).rstrip(".")
        ## Assign higher rank to unassigned rank (see (*))
        prev = lca_res[val_cols[0]]
        for rank in val_cols[1:]:
            if lca_res[rank]: prev = lca_res[rank]
            else: lca_res[rank] = "Unclassified."+prev
        lcadf[contig_id] = lca_res
    return pd.DataFrame.from_dict(lcadf,orient='index')

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Tabular file of contig_ids and taxids")
    parser.add_argument("-t", "--treefile", required=True,
            help="NCBI tree file (ncbi.tre)")
    parser.add_argument("-n", "--nodesfile", required=True,
            help="NCBI nodes file (nodes.dmp)")
    parser.add_argument("-m", "--mapfile", required=True,
            help="NCBI map file (ncbi.map)")
    
    args = parser.parse_args()
    
    ## Read the results from mash_parse.py
    result_df = pd.read_table(args.infile, sep='\t', index_col=0, header=0, names=['contig_id', 'taxid'])
    
    ## Get all unique taxids
    taxids = list(set(list(result_df['taxid'])))
    
    ## Load the NCBI files
    ncbi_tree = TreeNode.read(args.treefile)
    ncbi_map = pd.read_table(args.mapfile, index_col=0, header=None, names=['taxid','name'],usecols=[0,1])
    ranks_df = load_ranks(args.nodesfile)

    ## Limit assignments to these ranks
    val_cols = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    ## Assign ancestor taxids to each taxid found in the result file
    taxdf = assign_ancestors(taxids,ncbi_tree,ranks_df,val_cols)

    ## Assign LCAs for each contig
    lcadf = lca(result_df,taxdf,ncbi_map,val_cols)

    ## Write to file
    lcadf = lcadf.loc[:,val_cols]
    lcadf.to_csv(sys.stdout, header=True, sep="\t")
    
if __name__ == '__main__': main()
