#!/usr/bin/env python
###################################

from argparse import ArgumentParser
import sys,csv

def get_genome_sizes(f, genomes):
    hin = open(f, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    for row in hincsv: genomes[row[0]] = float(row[1])
    hin.close()
    return genomes

def read_blast(infile, i, length, sizenorm, mismatchid, besthit):
    samples = {}
    besthits = {}
    genomes = {}
    hin = open(infile, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    for row in hincsv:
        ## Define query, sample and subject
        query = row[0]
        subject = row[1]
        sample = query.split("=")[1]
        ## Remove '"' from sample id
        sample = sample.replace('"',"")
        genome = subject.split("=")[1]
        ## Check if sample is defined
        try: samples[sample]
        except KeyError: samples[sample] = {}
        ## Check identity
        if mismatchid: id = (1-float(row[4])/float(row[3]))*100
        else: id = float(row[2])
        if id<i: continue
        ## Check length
        l = int(row[3])
        if l < length: continue
        ## If counting best hits only, check if query has been counted already
        if besthit:
            try: 
                besthits[query]
                continue
            except KeyError: besthits[query] = ""
        ### Add counts for genome in sample
        try: samples[sample][genome]+=1
        except KeyError: samples[sample][genome] = 1
        ## Add genome to keys
        if not genome in genomes.keys(): genomes[genome] = 1
    hin.close()
    return samples,genomes

def write(outfile, samples, genomes):
    ## Write results
    if outfile: hout = open(outfile, 'w')
    else: hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    houtcsv.writerow(["Genome"]+sorted(samples.keys()))
    
    for genome in genomes.keys():
        out = [genome]
        size = genomes[genome]
        for sample in sorted(samples.keys()):
            try:
                hits = samples[sample][genome]
                if size>0: hits = float(hits)/float(size)
                out.append(hits)
            except KeyError: out.append(0.0)
        houtcsv.writerow(out)
    hout.close()

def main():
    parser = ArgumentParser(description='''This script produces tabular output from blast files for easier plotting with R.
    Sample ids should be attached to each query in the form of "sample_id=<sampleId>". 
    Genome ids should be attached to each subject in the form of "genome_id=<genomeId>".''')
    parser.add_argument("-i", "--infile", type=str, required=True, help="BLAST results file")
    parser.add_argument("-o", "--outfile", type=str, help="Table output")
    parser.add_argument("-I", "--id", type=float, default=99.0, help="Identity cutoff. Defaults to 99.0")
    parser.add_argument("-l", "--length", type=int, default=200, help="Alignment length cutoff. Defaults to 200 bp")
    parser.add_argument("-b", "--besthit", action="store_true", help="Count only best hits for each read")
    parser.add_argument("-s", "--sizenorm", type=str, help="Normalize number of hits by genome size. Specify TSV table with genomeids and genome sizes")
    parser.add_argument("-m", "--mismatchid", action="store_true", help="Calculate identity as number_of_mismatches/alignment_length (i.e. column 5/column 4")
    args = parser.parse_args()

    ## Read BLAST results
    samples, genomes = read_blast(args.infile, args.id, args.length, args.sizenorm, args.mismatchid, args.besthit)

    if args.sizenorm: genomes = get_genome_sizes(args.sizenorm, genomes)
    
    write(args.outfile, samples, genomes)

if __name__ == '__main__':
    main()
