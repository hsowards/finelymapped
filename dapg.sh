#!/bin/bash
module load plink
module load R
module load samtools

#arguments
#ARGUMENT=$@
rsid=$1
chr=$2
pos=$3
dist=$4
prefix=$5
causals=$6

#setting region
start=$(($pos-$dist))
end=$(($pos+$dist))

#test to make sure arguments were read properly, will display in slurm output
echo "rsid: $rsid"
echo "location: $chr:$start-$end"
echo "running DAP-G with $causals causals"

echo "Extracting 1000G LD Matrix for region"
#1000G
#LD
cd /data/Brown_lab/dap/LD_1000G #will be output to LD_1000G
bash /data/Brown_lab/vQTL/colocalization/LD_vQTL.sh $rsid $chr $pos $dist EUR $prefix /data/Brown_lab/vQTL/vQTL.env 

echo "Extracting summary statistics data for 1000G"
cd /data/Brown_lab/dap/sumstats
#$tabix ${vQTLfolder}/Melanocytes/Data/all_meta_gwas_2019.txt.gz ${Chr}:${Start}-${End} |awk -F "\t" -v OFS="\t" 'BEGIN{print "chr\tpos\tref\talt\trsnum\tpvalue\tzscore\teffect\tse"}{print $2,$3,$4,$5,$1,$6,$7,$8,$9}'  >${prefix}.GWAS.txt
tabix /data/Brown_lab/vQTL/Melanocytes/Data/all_meta_gwas_2019.txt.gz ${chr}:${start}-${end} |awk -F "\t" -v OFS="\t" 'BEGIN{print "chr\tpos\tref\talt\trsnum\tpvalue\tzscore\teffect\tse\tor\ta1\ta2\tn\tcase_n\tcontrol_n\tsample_n"}{print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}'  >${prefix}.GWAS.txt

echo "Extracting UKBB LD Matrix for region"
#UKBB
cd /data/Brown_lab/dap/LD_UKBB
plink --bfile /data/Brown_lab/vQTL/LD_reference/UKBB/UKBB_LD_bestguess_ref_panel_maf1e3_rsq_point3_HRC_SNPs_only --chr $chr --from-bp $start --to-bp $end --make-bed --out $prefix 

#Rscripts for data sets
echo "Aligning UKBB datasets for input into DAP-G"
Rscript /data/Brown_lab/dap/scripts/finemap_data_ukbb.R $chr $pos $dist $prefix
echo "Aligning 1000G datasets for input into DAP-G"
Rscript /data/Brown_lab/dap/scripts/finemap_data_1000G.R $chr $pos $dist $prefix

echo "Running DAP-G on the region using the UKBB LD reference panel"
cd /data/Brown_lab/dap/results
mkdir $prefix
cd $prefix
#dap-g, msize is num of causals
#UKBB
/home/sowardsha/modulefiles/dap/dap_src/dap-g -d_z /data/Brown_lab/dap/dap_input/sumstats/${prefix}_sumstats_UKBB.txt -d_ld /data/Brown_lab/dap/dap_input/ld/${prefix}_LD_UKBB.txt -t 2 -msize $causals > ${prefix}_${causals}_results_UKBB.txt

echo "Running DAP-G on the region using the 1000G LD reference panel"
#1000G
/home/sowardsha/modulefiles/dap/dap_src/dap-g -d_z /data/Brown_lab/dap/dap_input/sumstats/${prefix}_sumstats_1000G.txt -d_ld /data/Brown_lab/dap/dap_input/ld/${prefix}_LD_1000G.txt -t 2 -msize $causals > ${prefix}_${causals}_results_1000G.txt

