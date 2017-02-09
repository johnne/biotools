#!/usr/bin/env python
###################################

from ete3 import NCBITaxa
import pandas as pd
def main():
    ncbi = NCBITaxa()
    tree = ncbi.get_descendant_taxa(1,return_tree=True,intermediate_nodes=True)
    tree.write(outfile="ncbi.tre",format=8)
    
    r = tree.get_tree_root()
    taxids = [child.name for child in tree.get_descendants()]
    taxids.append(r.name)
    
    ncbi_map = pd.DataFrame.from_dict(ncbi.get_taxid_translator(taxids),orient='index')
    ncbi_map.to_csv("ncbi.map", sep="\t",header=False)

if __name__ == '__main__': main()
