#!/bin/bash

#input is a tab delimited arguments file with no header and the columns: rsid chr pos dist prefix causals
input="/data/Brown_lab/dap/Alkes.txt"

while IFS=$'\t' read -r line
do
echo $line
 sh /data/Brown_lab/dap/scripts/npztomatrix.sh $line
done < "$input"
