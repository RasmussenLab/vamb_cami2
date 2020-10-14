#!/usr/bin/env python

PPN = config.get("ppn", "")
MM_MEM = config.get("minimap_mem", "")
SAMPLES = config.get("sample_file", "")
SAMPLE_DATA = config.get("sample_data", "")
CONTIGS = config.get("contigs", "")
PROJECT = config.get("project", "")
VAMB_PARAMS = config.get("vamb_params", "")
TMP_DIR = config.get("tmpdir", "")
MIN_FASTA = config.get("min_fasta", "2000")
FLTFASTA = "fasta_filterlength.py"
NPZ = "ab_to_npz.py"

# read sample ids
fh_in = open(SAMPLES, 'r')
IDS = []
for line in fh_in:
    line = line.rstrip()
    IDS.append(line)
fh_in.close()

# read in sample2path
sample2path = {}
fh_in = open(SAMPLE_DATA, 'r')
for line in fh_in:
    line = line.rstrip()
    fields = line.split('\t')
    # only take samples in the subset that we are working on
    if fields[0] in IDS:
        sample2path[fields[0]] = fields[1]

rule all:
    input:
        expand("jgi/{sample}.cut.jgi", sample=IDS),
        "jgi_matrix/jgi.abundance.npz",
        "vamb/clusters.tsv",
        "metabat2/clusters.metabat.tsv"

rule filter_contigs:
    input:
        contigs = CONTIGS
    output:
        "contigs.flt.fna.gz"
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb", project=PROJECT
    threads:
        int(1)
    log:
        "log/filter_contigs/filter_contigs.log"
    shell:
        "{FLTFASTA} --i {input} --min {MIN_FASTA} | gzip -c > {output}"

rule index:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        mmi = "contigs.flt.mmi"
    params:
        walltime="864000", nodes="1", ppn="1", mem="90gb", project=PROJECT
    threads:
        int(1)
    log:
        "log/index/index.log"
    shell:
        "module load minimap2/2.17r941 samtools/1.10;"
        "minimap2 -I 12G -d {output} {input} 2> {log}"

rule dict:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        dict = "contigs.flt.dict"
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb", project=PROJECT
    threads:
        int(1)
    log:
        "log/dict/dict.log"
    shell:
        "module load samtools/1.10;"
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap:
    input:
        fq = lambda wildcards: sample2path[wildcards.sample],
        mmi = "contigs.flt.mmi",
        dict = "contigs.flt.dict"
    output:
        bam = temp("mapped/{sample}.bam")
    params:
        walltime="864000", nodes="1", ppn=PPN, mem=MM_MEM, project=PROJECT,
    threads:
        int(PPN)
    log:
        "log/map/{sample}.minimap.log"
    shell:
        "module load minimap2/2.17r941 samtools/1.10;"
        '''minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}'''

rule jgi:
    input:
        bam = "mapped/{sample}.bam"
        #aggregate_input
    output:
        jgi = temp("jgi/{sample}.raw.jgi")
    params:
        walltime="864000", nodes="1", ppn="1", mem="10gb", project=PROJECT
    threads:
        int(1)
    log:
        "log/jgi/{sample}.jgi"
    shell:
        "module unload R/3.6.1 gcc/8.2.0; module load perl/5.24.0 metabat/2.10.2;"
        "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input} 2>{log}"

rule cut_column1to3: 
    input:
        "jgi/%s.raw.jgi" % IDS[0] 
    output:
        "jgi/jgi.column1to3"
    params:
        walltime="86400", nodes="1", ppn="1", mem="1gb", project=PROJECT
    log:
        "log/jgi/column1to3"
    shell: 
        "cut -f1-3 {input} > {output} 2> {log}"

rule cut_column4to5:
    input:
        "jgi/{sample}.raw.jgi"
    output:
        "jgi/{sample}.cut.jgi"
    params:
        walltime="86400", nodes="1", ppn="1", mem="1gb", project=PROJECT
    log:
        "log/jgi_process/{sample}.cut.log"
    shell: 
        "cut -f1-3 --complement {input} > {output} 2> {log}"

rule paste_abundances:
    input:
        column1to3="jgi/jgi.column1to3",
        data=expand("jgi/{sample}.cut.jgi", sample=IDS)
    output:
        temp("jgi_matrix/jgi.abundance.dat") 
    params:
        walltime="86400", nodes="1", ppn="1", mem="1gb", project=PROJECT
    log:
        "log/jgi_process/paste_abundances.log"
    shell: 
        "paste {input.column1to3} {input.data} > {output} 2> {log}" 

rule ab_npz:
    input:
        "jgi_matrix/jgi.abundance.dat"
    output:
        "jgi_matrix/jgi.abundance.npz"
    params:
        walltime="86400", nodes="1", ppn="1", mem="100gb", project=PROJECT
    log:
        "log/jgi_process/ab2npz.log"
    shell: 
        "{NPZ} --i {input} --o {output} 2> {log}"  

rule vamb:
    input:
        jgi = "jgi_matrix/jgi.abundance.npz",
        contigs = "contigs.flt.fna.gz"
    output:
        "vamb/clusters.tsv",
        "vamb/latent.npz",
        "vamb/lengths.npz",
        "vamb/log.txt",
        "vamb/model.pt",
        "vamb/mask.npz",
        "vamb/tnf.npz"
    params:
        walltime="86400", nodes="1", ppn="4:gpus=1", mem="20gb", project=PROJECT
    log:
        "log/vamb/vamb.log"
    threads:
        int(4)
    shell:
        "rm -rf vamb;"
        "module load cuda/toolkit/10.1/10.1.243;"
        "vamb --outdir vamb --fasta {input.contigs} --rpkm {input.jgi} {VAMB_PARAMS}"



