## CoassemblyPipeline

Snakemake pipeline to generate and QC genome coassemblies from G&T-Seq


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

Example command to launch snakemake:

```
snakemake -p --config input=input.csv -j 20 --retries 1 --latency-wait 60 --snakefile CoassemblyPipeline.smk --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -c {cluster.c} --mem={cluster.memory} --job-name={cluster.J} --time={cluster.time} --exclude={cluster.exclude} --constraint=intel -o slurm_logs/slurm.{cluster.J}.%j.out"
```