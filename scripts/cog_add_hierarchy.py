#!/usr/bin/env python
###################################

from argparse import ArgumentParser
import sys, urllib2, csv

def write(cogs, infile, outfile):
    if outfile: hout = open(outfile, 'w')
    else: hout = sys.stdout
    hin = open(args.infile, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    houtcsv = csv.writer(hout, delimiter = '\t')
    for i,row in enumerate(hincsv):
        if i==0: 
            header= row[1:]
            out = ["Cat1","Cat2","COG"]+header
            houtcsv.writerow(out)
            continue
        cog = row[0]
        try: out = [cogs[cog]["firstcat"],cogs[cog]["secondcat"],cog+" "+cogs[cog]["name"]]+row[1:]
        except KeyError: out = [cog,cog,cog]+row[1:]
        houtcsv.writerow(out)
    hout.close()
    hin.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
            help="COG table for which to add hierarchy")
    #parser.add_argument("-H", "--hier", type=str,
    #        help="Supply a COG functional hierarchy flatfile. If none is supplied the script attempts to download one from ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt")
    parser.add_argument("-o", "--outfile", type=str,
            help="Write tab-delimited pathway hierarchy to outfile. Defaults to stdout")
    args = parser.parse_args()
    
    if not args.infile: sys.exit(parser.print_help())

    cogs = ParseWHOG()
    cogs = ParseFUN(cogs)
    write(cogs, args.infile, args.outfile)

def ParseWHOG():
    cogs = {}
    whogurl = "ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog"
    data = urllib2.urlopen(whogurl)
    data = data.read()
    data = data.split("\n")
    for line in data:
        if line=="": continue
        if not line[0]=="[": continue
        items = line.rsplit()
        cats = items[0]
        cog = items[1]
        name = " ".join(items[2:])
        cats = cats.replace("[","")
        cats = cats.replace("]","")
        cogs[cog] = {"name": name, "catletter": cats[0], "firstcat": "", "secondcat": ""}
    return cogs

def ParseFUN(cogs):
    f = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt")
    data = f.read()
    data = data.split("\n")
    cats = {}
    for line in data:
        line = line.lstrip()
        if line=="": continue
        if line[0]=="[":
            catletter = line[1]
            name = line.split("]")[1].lstrip()
            name = name.rstrip()
            cats[catletter] = {"secondcat": name, "firstcat":parent}
        else: parent = line
    for cog in cogs.keys():
        catletter = cogs[cog]["catletter"]
        firstcat = cats[catletter]["firstcat"]
        secondcat = cats[catletter]["secondcat"]
        cogs[cog]["firstcat"] = firstcat
        cogs[cog]["secondcat"] = secondcat
    return cogs

if __name__ == '__main__':
    main()
