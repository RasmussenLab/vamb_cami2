# GitHub repository for binning CAMI2 challenge dataset using VAMB

This is a snakemake file for running VAMB on the CAMI2 challenge datasets.

**IMPORTANT: Because CAMI datasets are pooled assemblies - this is not using the *multi-split* approach. If you are running your own data you *should* use the multi-split approach instead. I will put up a link here to a snakemake that does exactly that**

Edit the configuration file as necessary, the fields are:

```
configfile: the config-file itself
contigs: input contigs
minimap_mem: memory for minimap2
ppn: number of cores for minimap2
project: if running using qsub for accounting (-A option in qsub)
sample_file: text file with sample names to include (one per line)
sample_data: tab separated text file with sample_name _tab_ path/to/reads (it assumes interleaved reads)
tmpdir: temporary folder to write tmpstuff
vamb_params: extra parameters for vamb
```

### Example running
```
snakemake --cores 50 --configfile config.json --latency-wait 60
```

### Example running using qsub
```
snakemake --jobs 50 --configfile config.json --cluster "qsub -l walltime={params.walltime} -l nodes=1:ppn={params.ppn} -l mem={params.mem} -A cpr_10006 -W group_list=cpr_10006" --latency-wait 60
```
