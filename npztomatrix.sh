#!/bin/bash
prefix=$1

module load python
python3 /data/Brown_lab/dap/scripts/npztomatrix.py /data/Brown_lab/UKBB_LD_Alkes/${prefix}.npz /data/Brown_lab/UKBB_LD_Alkes/${prefix}.txt.gz
