GitHub repository for binning CAMI2 challenge dataset using VAMB

# example running
snakemake --cores 50 --configfile config.json --latency-wait 60

# example running using qsub
snakemake --jobs 50 --configfile config.json --cluster "qsub -l walltime={params.walltime} -l nodes=1:ppn={params.ppn} -l mem={params.mem} -A cpr_10006 -W group_list=cpr_10006" --latency-wait 60
