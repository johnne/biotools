#!/usr/bin/env python
###################################

from argparse import ArgumentParser
import sys,csv
def readData(infile):
    ## Read the hierarchy flatfile
    data = []
    if infile: 
        hin = open(args.infile, 'r')
        data = hin.read()
    else:
        import urllib2
        f = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt")
        data = f.read()
    data = data.split("\n")
    return data

def write(outfile, cats):
    if outfile: hout = open(args.outfile, 'w')
    else: hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    for cat, d in cats.iteritems(): houtcsv.writerow([cat,d["name"],d["hier"]])
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
            help="Supply a COG functional hierarchy flatfile. If none is supplied the script attempts to download one from ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt")
    parser.add_argument("-o", "--outfile", type=str,
            help="Write tab-delimited pathway hierarchy to outfile. Defaults to stdout")
    parser.add_argument("-I", "--include", type=str,
            help="Include only these pathways in the lookup")
    parser.add_argument("-E", "--exclude", type=str,
            help="Exclude these pathways in the lookup")
    args = parser.parse_args()
    exclude = []
    if args.exclude: exclude = ReadFile(args.exclude)
    if args.include: include = ReadFile(args.include)
    else: include = False

    data = readData(args.infile)
    cats = Parse(data, include, exclude)
    
    write(args.outfile, cats)

def ReadFile(f):
    d = {}
    hin = open(f, 'r')
    for line in hin: d[line.rstrip()] = ""
    hin.close()
    return d.keys()

def Parse(data, include, exclude):
    cats = {}
    for line in data:
        line = line.lstrip()
        if line=="": continue
        if line[0]=="[":
            cat = line[1]
            name = line.split("]")[1].lstrip()
            cats[cat] = {"name": name, "hier":parent}
        else: parent = line
    return cats

if __name__ == '__main__':
    main()
