#!/usr/bin/env python
###################################
import gzip as gz, csv, sys

def GetFilelist(i):
    import os
    ## If user supplies an existing hist file, use it
    if os.path.exists(i) and not "/dev/" in i: return [i]
    ## Otherwise, read the files
    files = []
    hin = open(i, 'r')
    for line in hin: files.append(line.rstrip())
    hin.close()
    return files

def WriteCov(samples, contigs, houtcsv):
    cover = {}
    houtcsv.writerow(["Contig"]+sorted(samples.keys()))    
    for contig in contigs.keys():
        out = [contig]
        for sample in sorted(samples.keys()):
            try: c = samples[sample][contig]
            except KeyError: c = 0.0
            try: cover[sample].append(c)
            except KeyError: cover[sample] = [c]
            out.append(c)
        houtcsv.writerow(out)
    return cover

def WriteStat(cover, houtcsv):
    import numpy as np
    houtcsv.writerow(["Sample","MeanCoverage","MedianCoverage","MaxCoverage","MinCoverage"])
    for sample,l in cover.iteritems():
        median_cov = np.median(l)
        max_cov = np.max(l)
        min_cov = np.min(l)
        mean_cov = np.mean(l)
        houtcsv.writerow([sample,mean_cov,median_cov,max_cov,min_cov])

def CalcCov(files):
    samples = {}
    contigs = {}
    cover = {}
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    for f in files:
        if ".gz" in f: hin = gz.open(f)
        else: hin = open(f)
        sample=(f.split("/")[-1]).split(".")[0]
        samples[sample] = {}
        cover[sample] = []
        hincsv = csv.reader(hin, delimiter = '\t')
        contig = ""
        cov = 0.0
        for i,row in enumerate(hincsv):
            if contig != row[0]:
                if len(files)==1:
                    if i==0:houtcsv.writerow(["Contig",sample])
                    else: houtcsv.writerow([contig,cov])
                else: samples[sample][contig] = cov
                contig = row[0]
                cov = 0.0
            contigs[contig] = ""
            depth = float(row[-4])
            depth_f = float(row[-1])
            cov += depth*depth_f
        hin.close()
        ## Handle last contig in the file
        if len(files)==1: houtcsv.writerow([contig,cov])
        else: samples[sample][contig] = cov
    if len(files)==1: sys.exit()
    cover = WriteCov(samples, contigs, houtcsv)
    hout.close()
    hout = sys.stderr
    houtcsv = csv.writer(hout, delimiter = '\t')
    WriteStat(cover, houtcsv)
    hout.close()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", type=str, required=True,
            help="Infile histogram(s) from bedtools coverage")
    args = parser.parse_args()
    files = GetFilelist(args.infiles)
    CalcCov(files)

if __name__ == "__main__": main()
