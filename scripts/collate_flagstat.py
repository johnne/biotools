#!/usr/bin/env python

from argparse import ArgumentParser
import sys

def get_files(f):
    files = []
    hin = open(f)
    for line in hin: files.append(line.rstrip())
    hin.close()
    return files

def collate_stats(f):
    r = {}
    hin = open(f)
    for line in hin:
        line = line.rstrip()
        if "in total (" in line: r["total"] = int(line.split(" ")[0])
        if "mapped (" in line: r["mapped"] = int(line.split(" ")[0])
        if "properly paired (" in line: r["prop-paired"] = int(line.split(" ")[0])
    hin.close()
    return r

def write(stats):
    import csv
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    houtcsv.writerow(["Sample","Total","Mapped","Mapped%","Properly-paired","Properly-paired%"])
    for n in sorted(stats.keys()):
        mapped_p = round(float(stats[n]["mapped"])/stats[n]["total"]*100,2)
        prop_p = round(float(stats[n]["prop-paired"])/stats[n]["total"]*100,2)
        houtcsv.writerow([n,stats[n]["total"],stats[n]["mapped"],mapped_p,stats[n]["prop-paired"],prop_p])
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
            help="Flagstat output")
    parser.add_argument("--files", type=str,
            help="List of flagstat output files to collate from")
    parser.add_argument("--no_extractname", action="store_true",
            help="Don't attempt to extract sample name from file name (defaults to False)")
    
    args = parser.parse_args()

    if not args.infile and not args.files: sys.exit(parser.print_help())
    
    if args.infile: files = [args.infile]
    else: files = get_files(args.files)

    stats = {}
    for f in files:
        name = f.split("/")[-1]
        if not args.no_extractname: name = name.split(".")[0]
        stats[name] = collate_stats(f)
    
    write(stats)

if __name__ == '__main__': 
    main()
