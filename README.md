## CoassemblyPipeline

Snakemake pipeline to generate and QC genome coassemblies from G&T-Seq


**THIS PIPELINE DOES NOT TRIM READS. IT EXPECTS THE INPUT READS TO HAVE ALREADY BEEN ADAPTER TRIMMED**


The input to this pipeline is a comma separated file. The first line specifies the name of the coassembly, the following lines contain 4 comma separated values to specify the gDNA and cDNA libraries:

- sample name
- library type (gDNA or cDNA)
- path to forward reads
- path to reverse reads

**Sample names should not contain periods**


Example:
`input.csv`

```
COASSEMBLY_NAME
Sample1,gDNA,Sample1_gDNA_R1.fq.gz,Sample1_gDNA_R2.fq.gz
Sample1,cDNA,Sample1_cDNA_R1.fg.gz,Sample1_cDNA_R2.fq.gz
Sample2,gDNA,Sample2_gDNA_R1.fq.gz,Sample2_gDNA_R2.fq.gz
Sample2,cDNA,Sample2_cDNA_R1.fg.gz,Sample2_cDNA_R2.fq.gz
Sample3,gDNA,Sample3_gDNA_R1.fq.gz,Sample3_gDNA_R2.fq.gz
Sample4,cDNA,Sample4_cDNA_R1.fq.gz,Sample4_cDNA_R2.fq.gz

```

Real example:

```
Diatom_Coassembly
S0002,gDNA,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0002_OX5002404_A69456/trimmed_reads/S0002_OX5002404_A69456_gDNA_R1.trimmed.fastq.gz,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0002_OX5002404_A69456/trimmed_reads/S0002_OX5002404_A69456_gDNA_R2.trimmed.fastq.gz
S0041,gDNA,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0041_OX5002443_A69495/trimmed_reads/S0041_OX5002443_A69495_gDNA_R1.trimmed.fastq.gz,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0041_OX5002443_A69495/trimmed_reads/S0041_OX5002443_A69495_gDNA_R2.trimmed.fastq.gz
S0042,cDNA,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0042_OX5002444_A69496/trimmed_reads/S0042_OX5002444_A69496_cDNA_R1.trimmed.fastq.gz,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0042_OX5002444_A69496/trimmed_reads/S0042_OX5002444_A69496_cDNA_R2.trimmed.fastq.gz
S0044,gDNA,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0044_OX5002446_A69498/trimmed_reads/S0044_OX5002446_A69498_gDNA_R1.trimmed.fastq.gz,/ei/.project-scratch/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/CB-GENANNO-476_DToL_Protists/Analysis/GENANNO-536_PL0325_deep_seq/S0044_OX5002446_A69498/trimmed_reads/S0044_OX5002446_A69498_gDNA_R2.trimmed.fastq.gz
```
