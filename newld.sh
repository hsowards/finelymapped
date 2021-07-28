#!/bin/bash

module load python

prefix=$1

#python3 npztotxt.py /data/Brown_lab/UKBB_LD_Alkes/${prefix}.npz /data/Brown_lab/UKBB_LD_Alkes/${prefix}.txt.gz
python3 npztotxt_pd.py /data/Brown_lab/UKBB_LD_Alkes/${prefix}.npz
#Rscript npztotxt.R /data/Brown_lab/UKBB_LD_Alkes/${prefix}.npz
