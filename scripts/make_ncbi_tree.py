#!/usr/bin/env python
###################################

from ete3 import NCBITaxa
ncbi = NCBITaxa()
tree = ncbi.get_descendant_taxa(1,return_tree=True,intermediate_nodes=True)
tree.write(outfile="ncbi.tre",format=8)
