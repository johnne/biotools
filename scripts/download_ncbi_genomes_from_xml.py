#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from collections import defaultdict
from ftplib import FTP

def parse(data):
    genbank_link = taxid = bioacc = ""
    for l in data:
        l = l.lstrip()
        items = l.split(">")
        if len(items)<3: continue
        key = items[0].lstrip("<")
        value = items[1].split("<")[0]
        if key=="FtpPath_GenBank": genbank_link = value
        elif key=="Taxid": taxid = value
        elif key=="BioprojectAccn": bioacc = value
        if taxid and genbank_link: return [taxid,bioacc,genbank_link]
    return [taxid,bioacc,genbank_link]

def read_lines(fh):
    data = []
    result = defaultdict(list)
    for line in fh:
        line = line.rstrip()
        #if line=="": continue
        data.append(line)
        if line=="</DocumentSummary>": 
            [taxid,bioacc,link] = parse(data)
            result[taxid+"."+bioacc]=[link]
            data = []
    return result

def download(key,link,outdir):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    genome_dir = "/".join(link.split("/")[3:])
    basename_dir = genome_dir.split("/")[-1]
    ftp.cwd(genome_dir)
    filename = basename_dir+"_genomic.fna.gz"
    sys.stderr.write("Downloading "+filename+" from server. Storing as "+outdir+"/"+key+".fasta.gz\n")
    ftp.retrbinary('RETR '+filename, open(outdir+"/"+key+".fasta.gz", 'wb').write)
    ftp.quit()
    
def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
            help="XML (assembly_result.xml) downloaded from NCBI")
    parser.add_argument("-o", "--outdir", default=".",
            help="Store downloaded files in outdir. Defaults to current directory")
    args = parser.parse_args()

    ## Read the genome and assembly info from the xml file
    sys.stderr.write("Reading XML file\n")
    with open(args.input) as fh: result = read_lines(fh)  
    
    ## Download
    for key,link in result.items(): 
        #print(key,link[0],sep="\t")
        download(key,link[0],args.outdir)
    
if __name__ == '__main__':
    main()
