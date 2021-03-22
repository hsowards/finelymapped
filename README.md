# finelymapped
dapg.sh: this will prep the data to run dapg on one locus, requires following arguments:
$ sbatch scripts/dapg.sh rsid chr pos dist prefix causals

dapg_loop.sh: this will run dapg.sh on multiple loci, inputted in a tab-delimited file with the following column format and no header:
rsid	chr	pos	dist	prefix	causals

finemap_data_1000G.R: alignment of summary statistics and LD reference matrix, including formatting for input into dap-g

