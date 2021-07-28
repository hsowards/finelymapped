#!/bin/bash
prefix=$1

module load R
#Rscript txttomatrix.R
python3 npztomatrix.py /data/Brown_lab/UKBB_LD_Alkes/${prefix}.npz /data/Brown_lab/dap/${prefix}.txt.gz
