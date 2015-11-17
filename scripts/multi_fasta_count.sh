#!/bin/bash

## Usage example: grep -c '>' *.fasta | multi_fasta_count.sh

sum=0

while read entry;
do
    n=${entry#*:}
    let sum=$sum+$n
done
echo $sum
