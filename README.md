## CoassemblyPipeline

Snakemake pipeline to generate and QC genome (co)assemblies from (G&T-Seq).


**THIS PIPELINE DOES NOT TRIM READS. IT EXPECTS THE INPUT READS TO HAVE ALREADY BEEN ADAPTER TRIMMED**


The input to this pipeline is a comma separated file. The first line specifies the name of the coassembly, the following lines contain 4 comma separated values to specify the gDNA and cDNA libraries:

- sample name
- library type (gDNA or cDNA)
- path to forward reads
- path to reverse reads

**Sample names should not contain periods**


Example input file: `input.csv`

```
COASSEMBLY_NAME
Sample1,gDNA,Sample1_gDNA_R1.fq.gz,Sample1_gDNA_R2.fq.gz
Sample1,cDNA,Sample1_cDNA_R1.fg.gz,Sample1_cDNA_R2.fq.gz
Sample2,gDNA,Sample2_gDNA_R1.fq.gz,Sample2_gDNA_R2.fq.gz
Sample2,cDNA,Sample2_cDNA_R1.fg.gz,Sample2_cDNA_R2.fq.gz
Sample3,gDNA,Sample3_gDNA_R1.fq.gz,Sample3_gDNA_R2.fq.gz
Sample4,cDNA,Sample4_cDNA_R1.fq.gz,Sample4_cDNA_R2.fq.gz

```

Example config file: `config.yaml`

```
# File listing input reads, better setting this at command line with "--config input=XXX.csv"
# Should use absolute paths
input:

# Optional tools to run
run_checkm: true
run_busco: true

# SPAdes parameters
kmers: "21,33,55,77"
min_scaffold_length: 1000 # for bbtools reformat

# BUSCO parameters
busco_version: 3
busco_database: XXXXX
augustus_species:

# MetaBat2 parameters
metabat2_min_contig: 1500

# Tiara parameters
tiara_min_length: 1000

# EukRep parameters
eukrep_min_length: 1000

# Blobtools parameters
diamond_database: XXXXX
megablast_database: XXXXX
nodesdb: XXXXX
n_chunks: 4 # chunk the fasta file into N splits for the megablastn search; randomly to balance as input fasta will be sorted by sequence length

# CAT parameters
cat_database: XXXXX
taxonomy_dir: XXXXX
diamond_path: XXXXX

# pr2 database
pr2_database: XXXXX

# cleanup options
cleanup_spades: true
cleanup_cat: true
```

Example command to launch snakemake:

```
snakemake -p --config input=input.csv -j 20 --retries 1 --latency-wait 60 --snakefile CoassemblyPipeline.smk --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time} --exclude={cluster.exclude} -o slurm_logs/slurm.{cluster.J}.%j.out"
```

---

This pipeline runs:

- genome assembly using [SPAdes](https://github.com/ablab/spades)
- annotation of rRNA genes using [barrnap](https://github.com/tseemann/barrnap), followed by classification based on comparison with the [pr2](https://github.com/pr2database/pr2database) database
- taxonomic classification using:
	- [Blobtools](https://github.com/DRL/blobtools)
	- [CAT](https://github.com/dutilh/CAT)
	- [EukRep](https://github.com/patrickwest/EukRep)
	- [Tiara](https://github.com/ibe-uw/tiara/)
- assembly stats using [QUAST](https://github.com/ablab/quast)
- coverage stats based on mapping reads with [minimap2](https://github.com/lh3/minimap2) and [QualiMap](http://qualimap.conesalab.org/)
- contig binning using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat)
- [BUSCO](https://gitlab.com/ezlab/busco)
- [CheckM](https://github.com/Ecogenomics/CheckM)
- Optionally, aligns RNA-Seq reads using [hisat2](https://github.com/DaehwanKimLab/hisat2) and assembles transcripts using [StringTie2](https://github.com/gpertea/stringtie)
