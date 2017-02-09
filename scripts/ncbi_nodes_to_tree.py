#!/usr/bin/env python

import pandas as pd, sys, networkx as nx
from argparse import ArgumentParser
from skbio.tree import TreeNode, MissingNodeError

def create_graph(nodes_file):
    treegraph = nx.DiGraph()
    with open(nodes_file) as nodes_f:
        for line in nodes_f:
            line = line.strip()
            data_fields = line.split('\t')
            taxid = data_fields[0]
            parent = data_fields[2]
            treegraph.add_edge(parent,taxid)
    return treegraph

def main():
    parser = ArgumentParser()
    parser.add_argument("-n", "--nodesfile", required=True,
            help="NCBI nodes file (nodes.dmp)")

    args = parser.parse_args()

    treegraph = create_graph(args.nodesfile)

    
