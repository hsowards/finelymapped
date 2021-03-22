#!/bin/bash

#input is a tab delimited arguments file with no header and the columns: rsid chr pos dist prefix causals
input="/data/Brown_lab/dap/args.txt"

#this loop runs each row (loci) from the command file through dapg.sh
while IFS=$'\t' read -r line
do
 sh /data/Brown_lab/dap/scripts/dapg.sh $line
done < "$input"
