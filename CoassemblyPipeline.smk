# Snakemake pipeline to generate genome coassemblies from G&T-Seq data
# 
# Jamie McGowan, 2022
# <jamie.mcgowan@earlham.ac.uk>

import sys
from os import mkdir
from os.path import join, basename, isdir, normpath, realpath
from os import getcwd
from re import split

configfile: "config.yaml"

assembly_name = ""
gDNA_libraries = {}
cDNA_libraries = {}

# print("Coassembly pipeline")
# print()
# print("Current working directory:", getcwd())

with open(config["input"]) as f:
    line = f.readline().strip().strip(",")

    if "," in line or " " in line:
        print("Error line 1 in", config["input"], "should specify co-assembly name")
        sys.exit(0)
    
    assembly_name = line

    for line in f:
        if len(line.strip()) == 0:
            continue

        line = [l.strip() for l in line.strip().split(",")]

        if len(line) != 4:
            print("Error wrong intput file format to specify library in", config["input"])
            print(line)
            sys.exit(0)

        sample_name = line[0]
        sample_type = line[1]
        r1_path = realpath(line[2])
        r2_path = realpath(line[3])

        if sample_type.upper() == "GDNA":
            if sample_name not in gDNA_libraries:
                gDNA_libraries[sample_name] = [r1_path, r2_path]
            else:
                print("Error", sample_name, "gDNA", "has been specified multiple times")
                sys.exit()
        elif sample_type.upper() == "CDNA":
            if sample_name not in cDNA_libraries:
                cDNA_libraries[sample_name] = [r1_path, r2_path]
            else:
                print("Error", sample_name, "cDNA", "has been specified multiple times")
                sys.exit(0)
        else:
            print("Error, don't understand sample type:", sample_type)
            print("Accepted values are gDNA and cDNA")
            sys.exit(0)

# print()
# print("User has specified", len(gDNA_libraries), "gDNA libraries")
# for sample in gDNA_libraries:
    # print(sample, gDNA_libraries[sample])

# print()
# print("User has specified", len(cDNA_libraries), "cDNA libraries")
# for sample in cDNA_libraries:
    # print(sample, cDNA_libraries[sample])

if len(gDNA_libraries) == 0:
    print("Error. You must specify at least one gDNA library")
    sys.exit(0)

gDNA_samples = list(gDNA_libraries)
cDNA_samples = list(cDNA_libraries)

# print()
# print("Assembly name:", assembly_name)
# print("Min scaffold length:", config["min_scaffold_length"])

filtered_scaffolds_filename = join(assembly_name, "assembly", "scaffolds." + str(config["min_scaffold_length"]) + ".fasta")

TARGETS = []
# final summary
TARGETS.append(join(assembly_name, assembly_name + "_summary.tsv"))
#quast_report
TARGETS.append(join(assembly_name, "quast", "report.txt"))
#infoseq_summary
TARGETS.append(join(assembly_name, "assembly", "scaffolds." + str(config["min_scaffold_length"]) + ".tsv"))
#rRNA genes from barrnap
TARGETS.append(join(assembly_name, "rrna", assembly_name + ".gDNA.rrna.blast.top.tsv"))
#gDNA_depth
TARGETS.append(join(assembly_name, "gDNA_alignments", "depth.tsv"))
#hisat2_alignments
TARGETS.append(expand(join(assembly_name, "cDNA_alignments", "{sample}.cDNA.bam"), sample = cDNA_samples))
#metabat2_bins 
TARGETS.append(expand(join(assembly_name, "metabat2", "{coverage}", "bin_members.tsv"), coverage = ["coverage", "no_coverage"]))
#quast_report_bins
TARGETS.append(expand(join(assembly_name, "metabat2", "{coverage}", "quast", "report.txt"), coverage = ["coverage", "no_coverage"]))
#tiara_classification
TARGETS.append(join(assembly_name, "tiara", assembly_name + ".tiara.tsv"))
#eukrep_classification
TARGETS.append(join(assembly_name, "eukrep", assembly_name + ".eukrep.tsv"))
#checkm_done
TARGETS.append(expand(join(assembly_name, "checkm", "{coverage}", "checkm.done"), coverage = ["coverage", "no_coverage"]))
#blobtools_classification
TARGETS.append(join(assembly_name, "blobtools", assembly_name + ".blobDB.bestsum.table.txt"))
#CAT_summary
TARGETS.append(join(assembly_name, "CAT", assembly_name + ".CAT.summary.txt"))

#busco_done
if config["run_busco"] == True:
    TARGETS.append(join(assembly_name, "assembly", "busco.done"))

#checkm_done
if config["run_checkm"] == True:
    TARGETS.append(expand(join(assembly_name, "checkm", "{coverage}", "checkm.done"), coverage = ["coverage", "no_coverage"]))

if len(gDNA_samples) > 1:
    #qualimap_comparison_report
    TARGETS.append(join(assembly_name, "gDNA_alignments", "qualimap", "comparison", "multisampleBamQcReport.html"))
    #idxstats_gDNA
    TARGETS.append(expand(join(assembly_name, "gDNA_alignments", "{sample}.gDNA.idxstats.tsv"), sample = gDNA_samples + ["merged"]))
    #qualimap_reports
    TARGETS.append(expand(join(assembly_name, "gDNA_alignments", "qualimap", "{sample}", "qualimapReport.html"), sample = gDNA_samples + ["merged"]))
else:
    #idxstats_gDNA
    TARGETS.append(expand(join(assembly_name, "gDNA_alignments", "{sample}.gDNA.idxstats.tsv"), sample = gDNA_samples))
    #qualimap_reports
    TARGETS.append(expand(join(assembly_name, "gDNA_alignments", "qualimap", "{sample}", "qualimapReport.html"), sample = gDNA_samples))

if len(cDNA_samples) > 1:
    #merged_hisat2_alignment
    TARGETS.append(join(assembly_name, "cDNA_alignments", "merged.cDNA.bam"))
    #idxstats_cDNA
    TARGETS.append(expand(join(assembly_name, "cDNA_alignments", "{sample}.cDNA.idxstats.tsv"), sample = cDNA_samples + ["merged"]))
    #stringtie_assembie
    TARGETS.append(expand(join(assembly_name, "cDNA_alignments", "{sample}.stringtie.gtf"), sample = cDNA_samples + ["merged"]))
else:
    #idxstats_cDNA
    TARGETS.append(expand(join(assembly_name, "cDNA_alignments", "{sample}.cDNA.idxstats.tsv"), sample = cDNA_samples))
    #stringtie_assembie
    TARGETS.append(expand(join(assembly_name, "cDNA_alignments", "{sample}.stringtie.gtf"), sample = cDNA_samples))

rule all:
    input: TARGETS

# Writes a yaml file for SPAdes to speficy file paths to gDNA libraries
rule create_yaml_file:
    output:
        yaml_file = join(assembly_name, "input_dataset.yaml")
    run:
        with open(output.yaml_file, "w") as fo:
            
            fo.write('[\n')
            fo.write('   {\n')
            fo.write('      orientation: "fr",\n')
            fo.write('      type: "paired-end",\n')
            fo.write('      right reads: [\n')
            
            n = 0
            for lib in gDNA_libraries:
                n+= 1
                
                if n == len(gDNA_libraries):
                    fo.write('         "' + gDNA_libraries[lib][1] + '"\n')
                else:
                    fo.write('         "' + gDNA_libraries[lib][1] + '",\n')
            fo.write('      ],\n')
            fo.write('      left reads: [\n')

            n = 0
            for lib in gDNA_libraries:
                n += 1
                if n == len(gDNA_libraries):
                    fo.write('         "' + gDNA_libraries[lib][0] + '"\n')
                else:
                    fo.write('         "' + gDNA_libraries[lib][0] + '",\n')
            fo.write('      ]\n')
            fo.write('   }\n')
            fo.write(']\n')

rule spades:
    input:
        yaml_file = join(assembly_name, "input_dataset.yaml")
    params:
        kmers = config["kmers"],
        assembly_name = assembly_name,
        output_directory = join(assembly_name, "assembly")
    output:
        scaffolds = join(assembly_name, "assembly", "scaffolds.fasta")
    threads: 8
    log: join(assembly_name, "logs", "spades.log")
    benchmark: join(assembly_name, "benchmarks", "spades.tsv")
    shell: "/ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/SPAdes-3.15.5-Linux/bin/spades.py --dataset {input.yaml_file} --sc -k {params.kmers} --threads {threads} -o {params.output_directory}"

rule filter_assembly_length:
    input:
        scaffolds = join(assembly_name, "assembly", "scaffolds.fasta")
    output:
        filtered_scaffolds = filtered_scaffolds_filename
    params:
        min_scaffold_length = config["min_scaffold_length"]
    threads: 1
    log: join(assembly_name, "logs", "reformat.log")
    benchmark: join(assembly_name, "benchmarks", "reformat.tsv")
    shell: "reformat.sh in={input.scaffolds} out={output.filtered_scaffolds} minlength={params.min_scaffold_length}"

rule quast:
    input:
        scaffolds = join(assembly_name, "assembly", "scaffolds.fasta"),
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        quast_report = join(assembly_name, "quast", "report.txt")
    params:
        output_directory = directory(join(assembly_name, "quast"))
    threads: 1
    log: join(assembly_name, "logs", "quast.log")
    benchmark: join(assembly_name, "benchmarks", "quast.tsv")
    shell: "quast -t {threads} -o {params.output_directory} {input}  > {log} 2>&1"

rule infoseq:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        infoseq_summary = join(assembly_name, "assembly", "scaffolds." + str(config["min_scaffold_length"]) + ".tsv")
    threads: 1
    log: join(assembly_name, "logs", "infoseq.log")
    benchmark: join(assembly_name, "benchmarks", "infoseq.tsv")
    shell: "infoseq -sequence {input.filtered_scaffolds} -auto -nocolumns -delimiter '\t' -only -name -length -pgc -outfile {output.infoseq_summary}"

rule minimap2:
    input:
        filtered_scaffolds = filtered_scaffolds_filename,
        r1 = lambda wildcards: gDNA_libraries[wildcards.sample][0],
        r2 = lambda wildcards: gDNA_libraries[wildcards.sample][1],
    output:
        alignment = join(assembly_name, "gDNA_alignments", "{sample}.gDNA.bam")
    threads: 4
    log: join(assembly_name, "logs", "minimap2_{sample}.log")
    benchmark: join(assembly_name, "benchmarks", "minimap2_{sample}.tsv")
    shell: "minimap2 -t {threads} -ax sr {input.filtered_scaffolds} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output.alignment} && samtools index -@ {threads} {output.alignment}"

rule merge_minimap2:
    input:
        alignments = expand(join(assembly_name, "gDNA_alignments", "{sample}.gDNA.bam"), sample = gDNA_samples)
    output:
        merged_alignment = join(assembly_name, "gDNA_alignments", "merged.gDNA.bam"),
        index = join(assembly_name, "gDNA_alignments", "merged.gDNA.bam.bai")
    threads: 2
    log: join(assembly_name, "logs", "samtools_merge_minimap2.log")
    benchmark: join(assembly_name, "benchmarks", "samtools_merge_minimap2.tsv")
    shell: "samtools merge -@ {threads} {output.merged_alignment} {input.alignments} && samtools index -@ {threads} {output.merged_alignment}"

rule idxstats_minimap2:
    input:
        alignment = join(assembly_name, "gDNA_alignments", "{sample}.gDNA.bam")
    output:
        idxstats = join(assembly_name, "gDNA_alignments", "{sample}.gDNA.idxstats.tsv")
    threads: 1
    log: join(assembly_name, "logs", "samtools_idxstats_gDNA_{sample}.log")
    benchmark: join(assembly_name, "benchmarks", "samtools_idxstats_gDNA_{sample}.tsv")
    shell: "samtools idxstats -@ {threads} {input.alignment} > {output.idxstats}"

rule qualimap_bamqc: 
    input:
        alignment = join(assembly_name, "gDNA_alignments", "{sample}.gDNA.bam")
    output:
        report = join(assembly_name, "gDNA_alignments", "qualimap", "{sample}", "qualimapReport.html")
    params:
        output_directory = join(assembly_name, "gDNA_alignments", "qualimap", "{sample}")
    threads: 2
    shell: "qualimap bamqc -bam {input.alignment} -outdir {params.output_directory} -outformat PDF:HTML -nt {threads} --java-mem-size=20G"

rule qualimap_multi_bamqc_input:
    input:
        reports = expand(join(assembly_name, "gDNA_alignments", "qualimap", "{sample}", "qualimapReport.html"), sample = gDNA_samples)
    output:
        qualimap_multi_bamqc_input = join(assembly_name, "gDNA_alignments", "qualimap", "report_paths.tsv")
    params:
        samples = gDNA_samples,
        reports_directories = expand(join(assembly_name, "gDNA_alignments", "qualimap", "{sample}"), sample = gDNA_samples)
    threads: 1
    run:
        fo = open(output.qualimap_multi_bamqc_input, "w")
        for sample, report in zip(params.samples, params.reports_directories):
            fo.write("\t".join([sample, realpath(report)]) + "\n")
        fo.close()

rule qualimap_multi_bamqc:
    input:
        qualimap_multi_bamqc_input = join(assembly_name, "gDNA_alignments", "qualimap", "report_paths.tsv")
    output:
        report = join(assembly_name, "gDNA_alignments", "qualimap", "comparison", "multisampleBamQcReport.html")
    params:
        output_directory = join(assembly_name, "gDNA_alignments", "qualimap", "comparison")
    threads: 2
    shell: "qualimap multi-bamqc -d {input.qualimap_multi_bamqc_input} -outdir {params.output_directory} -outformat PDF:HTML"

rule hisat2_build:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        hisat2_index_done = join(assembly_name, "cDNA_alignments", "hisat2_index", "hisat2_build.done")
    params:
        index_prefix = join(assembly_name, "cDNA_alignments", "hisat2_index", assembly_name)
    threads: 4
    log: join(assembly_name, "logs", "hisat2_build.log")
    benchmark: join(assembly_name, "benchmarks", "hisat2_build.tsv")
    shell: "hisat2-build -p {threads} {input.filtered_scaffolds} {params.index_prefix} && touch {output.hisat2_index_done}"

rule hisat2:
    input:
        hisat2_index_done = join(assembly_name, "cDNA_alignments", "hisat2_index", "hisat2_build.done"),
        r1 = lambda wildcards: cDNA_libraries[wildcards.sample][0],
        r2 = lambda wildcards: cDNA_libraries[wildcards.sample][1]
    params:
        index_prefix = join(assembly_name, "cDNA_alignments", "hisat2_index", assembly_name)
    output:
        alignment = join(assembly_name, "cDNA_alignments", "{sample}.cDNA.bam"),
        summary = join(assembly_name, "cDNA_alignments", "{sample}.hisat2_summary.txt")
    threads: 4
    log: join(assembly_name, "logs", "hisat2_{sample}.log")
    benchmark: join(assembly_name, "benchmarks", "hisat2_{sample}.tsv")
    shell: 'hisat2 -q -p {threads} --dta --summary-file {output.summary} -t -x {params.index_prefix} -1 {input.r1} -2 {input.r2} | samtools sort -@ {threads} -O BAM -o {output.alignment} && samtools index -@ {threads} {output.alignment}'

rule merge_hisat2:
    input:
        alignments = expand(join(assembly_name, "cDNA_alignments", "{sample}.cDNA.bam"), sample = cDNA_samples)
    output:
        merged_alignment = join(assembly_name, "cDNA_alignments", "merged.cDNA.bam")
    threads: 2
    log: join(assembly_name, "logs", "samtools_merge_hisat2.log")
    benchmark: join(assembly_name, "benchmarks", "samtools_merge_hisat2.tsv")
    shell: "samtools merge -@ {threads} {output.merged_alignment} {input.alignments} && samtools index -@ {threads} {output.merged_alignment}"

rule idxstats_hisat2:
    input:
        alignment = join(assembly_name, "cDNA_alignments", "{sample}.cDNA.bam")
    output:
        idxstats = join(assembly_name, "cDNA_alignments", "{sample}.cDNA.idxstats.tsv")
    threads: 1
    log: join(assembly_name, "logs", "samtools_idxstats_cDNA_{sample}.log")
    benchmark: join(assembly_name, "benchmarks", "samtools_idxstats_cDNA_{sample}.tsv")
    shell: "samtools idxstats -@ {threads} {input.alignment} > {output.idxstats}"

rule stringtie:
    input:
        alignment = join(assembly_name, "cDNA_alignments", "{sample}.cDNA.bam")
    output:
        stringtie_assembly = join(assembly_name, "cDNA_alignments", "{sample}.stringtie.gtf")
    threads: 1
    log: join(assembly_name, "logs", "stringtie_{sample}.log")
    benchmark: join(assembly_name, "benchmarks", "stringtie_{sample}.tsv")
    shell: "stringtie -o {output.stringtie_assembly} {input.alignment}"

rule jgi_summarize_bam_contig_depths:
    input:
        alignments = expand(join(assembly_name, "gDNA_alignments", "{sample}.gDNA.bam"), sample = gDNA_samples)
    output:
        depth = join(assembly_name, "gDNA_alignments", "depth.tsv")
    threads: 1
    log: join(assembly_name, "logs", "jgi_summarize_bam_contig_depths.log")
    benchmark: join(assembly_name, "benchmarks", "jgi_summarize_bam_contig_depths.tsv")
    shell: "/ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/metabat2-2.15/bin/jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.alignments}"

# Run metabat2 without coverage information (i.e. tetra-nucleotide frequencies only)
rule metabat2:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        bin_members  = join(assembly_name, "metabat2", "no_coverage", "bin")
    params:
        min_contig = config["metabat2_min_contig"],
        metabat2_output = join(assembly_name, "metabat2", "no_coverage", "bin")
    threads: 2
    log: join(assembly_name, "logs", "metabat2_no_coverage.log")
    benchmark: join(assembly_name, "benchmarks", "metabat2_no_coverage.tsv")
    shell: "/ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/metabat2-2.15/bin/metabat2 -i {input.filtered_scaffolds} -o {params.metabat2_output} -m {params.min_contig} --numThreads {threads} --saveCls --unbinned"

# Run metabat2 with coverage information
rule metabat2_coverage:
    input:
        filtered_scaffolds = filtered_scaffolds_filename,
        depth = join(assembly_name, "gDNA_alignments", "depth.tsv")
    output:
        bin_members  = join(assembly_name, "metabat2", "coverage", "bin")
    params:
        min_contig = config["metabat2_min_contig"],
        metabat2_output = join(assembly_name, "metabat2", "coverage", "bin")
    threads: 2
    log: join(assembly_name, "logs", "metabat2_coverage.log")
    benchmark: join(assembly_name, "benchmarks", "metabat2_coverage.tsv")
    shell: "/ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/metabat2-2.15/bin/metabat2 -a {input.depth} -i {input.filtered_scaffolds} -o {params.metabat2_output} -m {params.min_contig} --numThreads {threads} --saveCls --unbinned"

rule parse_metabat2:
    input: 
        bin_members = join(assembly_name, "metabat2", "{coverage}", "bin")
    output:
        bin_members_parsed  = join(assembly_name, "metabat2", "{coverage}", "bin_members.tsv")
    params:
        metabat2_output_directory = join(assembly_name, "metabat2", "{coverage}"),
        metabat2_output_prefix = 'bin'
    shell: 'parse_metabat2.py {params.metabat2_output_directory} {params.metabat2_output_prefix} > {output.bin_members_parsed}'

rule quast_bins:
    input:
        bin_members = join(assembly_name, "metabat2", "{coverage}", "bin")
    output:
        quast_report_bins = join(assembly_name, "metabat2", "{coverage}", "quast", "report.txt")
    params:
        bins = join(assembly_name, "metabat2", "{coverage}", "*.fa"),
        output_directory = join(assembly_name, "metabat2", "{coverage}", "quast")
    threads: 2
    log: join(assembly_name, "logs", "quast_bins_{coverage}.log")
    benchmark: join(assembly_name, "benchmarks", "quast_bins_{coverage}")
    shell: "quast -t {threads} -o {params.output_directory} {params.bins}  > {log} 2>&1"

rule tiara:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        tiara_classification = join(assembly_name, "tiara", assembly_name + ".tiara.tsv")
    params:
        min_length = config["tiara_min_length"],
        tf = 'all'
    threads: 2
    log: join(assembly_name, "logs", "tiara.log")
    benchmark: join(assembly_name, "benchmarks", "tiara")
    shell: 'source tiara-1.0.1 && tiara -i {input.filtered_scaffolds} -o {output.tiara_classification} -m {params.min_length} --tf {params.tf} -t {threads} -v > {log} 2>&1'

rule eukrep:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        eukrep_eukaryotic = join(assembly_name, "eukrep", assembly_name + ".eukaryotic.fasta"),
        eukrep_prokaryotic = join(assembly_name, "eukrep", assembly_name + ".prokaryotic.fasta")
    params:
        min_length = config["eukrep_min_length"]
    threads: 2
    log: join(assembly_name, "logs", "eukrep.log")
    benchmark: join(assembly_name, "benchmarks", "eukrep")
    shell: 'source eukrep-0.6.6 && EukRep -i {input.filtered_scaffolds} -o {output.eukrep_eukaryotic} --min {params.min_length} --prokarya {output.eukrep_prokaryotic}'

rule parse_eukrep:
    input:
        eukrep_eukaryotic = join(assembly_name, "eukrep", assembly_name + ".eukaryotic.fasta"),
        eukrep_prokaryotic = join(assembly_name, "eukrep", assembly_name + ".prokaryotic.fasta")
    output:
        eukrep_classification = join(assembly_name, "eukrep", assembly_name + ".eukrep.tsv")
    threads: 1
    shell: 'parse_eukrep.py {input.eukrep_eukaryotic} {input.eukrep_prokaryotic} > {output.eukrep_classification}'

rule checkm:
    input:
        bin_members  = join(assembly_name, "metabat2", "{coverage}", "bin")
    output:
        checkm_done = join(assembly_name, "checkm", "{coverage}", "checkm.done")
    params:
        bins = join(assembly_name, "metabat2", "{coverage}"),
        output_dir = join(assembly_name, "checkm", "{coverage}"),
        checkm_out = join(assembly_name, "checkm", "{coverage}", "checkm.out")
    threads: 8
    log: join(assembly_name, "logs", "checkm_{coverage}.log")
    benchmark: join(assembly_name, "benchmarks", "checkm_{coverage}")
    shell: 'checkm lineage_wf {params.bins} {params.output_dir} -x fa -t {threads} > {params.checkm_out} && touch {output.checkm_done}'

rule barrnap:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        rrna = join(assembly_name, "rrna", assembly_name + ".gDNA.rrna.fasta")
    params:
        kingdom = 'euk'
    threads: 2
    log: join(assembly_name, "logs", "barrnap.log")
    benchmark: join(assembly_name, "benchmarks", "barrnap.tsv")
    shell: 'source barrnap-0.9 && barrnap --threads {threads} --kingdom {params.kingdom} --outseq {output.rrna} {input.filtered_scaffolds}'

rule pr2_blast:
    input:
        rrna = join(assembly_name, "rrna", assembly_name + ".gDNA.rrna.fasta")
    output:
        hits = join(assembly_name, "rrna", assembly_name + ".gDNA.rrna.blast.tsv"),
        best_hits = join(assembly_name, "rrna", assembly_name + ".gDNA.rrna.blast.top.tsv")
    params:
        database = config["pr2_database"]
    threads: 2
    log: join(assembly_name, "logs", "pr2_blast.log")
    benchmark: join(assembly_name, "benchmarks", "pr2_blast.tsv")
    shell: 'blastn -query {input.rrna} -db {params.database} -outfmt 6 -out {output.hits} -num_threads {threads} && print_best_blast_hits.py {output.hits} > {output.best_hits}'

rule diamond:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        hits = join(assembly_name, "blobtools", assembly_name + ".diamond.out")
    params:
        db = config["diamond_database"],
        outfmt= '6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sskingdoms skingdoms sphylums sscinames'
    threads: 16
    log: join(assembly_name, "logs", "diamond.log")
    benchmark: join(assembly_name, "benchmarks", "diamond.tsv")
    shell: 'source diamond-2.0.14 && /usr/bin/time -v diamond blastx --query {input.filtered_scaffolds} --sensitive --max-target-seqs 1 --evalue 1e-25 --threads {threads} --outfmt {params.outfmt} --db {params.db} > {output.hits}'

rule megablast:
    input:
        filtered_scaffolds = filtered_scaffolds_filename
    output:
        hits = join(assembly_name, "blobtools", assembly_name + ".blastn.out")
    params:
        db = config["megablast_database"],
        outfmt = '6 qseqid staxids bitscore std'
    threads: 6
    log: join(assembly_name, "logs", "megablast.log")
    benchmark: join(assembly_name, "benchmarks", "megablast.tsv")
    shell: 'blastn -task megablast -db {params.db} -query {input.filtered_scaffolds} -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads {threads} -out {output.hits} -outfmt "{params.outfmt}"'

rule blobtools:
    input:
        filtered_scaffolds = filtered_scaffolds_filename,
        merged_bam = join(assembly_name, "gDNA_alignments", "merged.gDNA.bam"),
        index = join(assembly_name, "gDNA_alignments", "merged.gDNA.bam.bai"),
        diamond_hits = join(assembly_name, "blobtools", assembly_name + ".diamond.out"),
        blastn_hits = join(assembly_name, "blobtools", assembly_name + ".blastn.out")
    output: 
        blobdb = join(assembly_name, "blobtools", assembly_name + ".blobDB.json"),
        view = join(assembly_name, "blobtools", assembly_name + ".blobDB.bestsum.table.txt")
    params:
        prefix = join(assembly_name, "blobtools", assembly_name),
        prefix_directory = join(assembly_name, "blobtools/"),
        nodesdb = config["nodesdb"]
    threads: 1
    log: join(assembly_name, "logs", "blobtools.log")
    benchmark: join(assembly_name, "benchmarks", "blobtools.tsv")
    shell:
        """
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools create -i {input.filtered_scaffolds} -b {input.merged_bam} -t {input.diamond_hits} -t {input.blastn_hits} -o {params.prefix} --db {params.nodesdb} > {log} 2>&1
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools view -i {output.blobdb} -o {params.prefix_directory} -r all >> {log} 2>&1
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools plot -i {output.blobdb} -o {params.prefix_directory} -r superkingdom >> {log} 2>&1
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools plot -i {output.blobdb} -o {params.prefix_directory} -r phylum >> {log} 2>&1
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools plot -i {output.blobdb} -o {params.prefix_directory} -r order >> {log} 2>&1
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools plot -i {output.blobdb} -o {params.prefix_directory} -r family >> {log} 2>&1
        singularity exec /hpc-home/mcgowan/software/singularity_testing/blobtools.img blobtools plot -i {output.blobdb} -o {params.prefix_directory} -r genus >> {log} 2>&1
        """

rule CAT_classify:
    input:
        filtered_scaffolds = filtered_scaffolds_filename,
    output:
        CAT_classification = join(assembly_name, "CAT", assembly_name + ".CAT.contig2classification.officialnames.txt"),
        CAT_summary = join(assembly_name, "CAT", assembly_name + ".CAT.summary.txt")
    params:
        output_prefix = assembly_name + ".CAT",
        output_directory = join(assembly_name, "CAT"),
        CAT_database = config["cat_database"],
        taxonomy_dir = config["taxonomy_dir"],
        diamond_path = config["diamond_path"]
    threads: 8
    log: join(assembly_name, "logs", "CAT.log")
    benchmark: join(assembly_name, "benchmarks", "CAT.tsv")
    shell: 
        """
        singularity exec /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/CAT-5.2/CAT.img CAT contigs -o {params.output_prefix} -c {input.filtered_scaffolds} -d {params.CAT_database} -t {params.taxonomy_dir} -n {threads} --path_to_diamond {params.diamond_path} --I_know_what_Im_doing --top 25
        singularity exec /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/CAT-5.2/CAT.img CAT add_names -i {params.output_prefix}.contig2classification.txt -o {params.output_prefix}.contig2classification.names.txt -t {params.taxonomy_dir}
        singularity exec /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/CAT-5.2/CAT.img CAT add_names -i {params.output_prefix}.contig2classification.txt -o {params.output_prefix}.contig2classification.officialnames.txt -t {params.taxonomy_dir} --only_official
        singularity exec /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/CAT-5.2/CAT.img CAT add_names -i {params.output_prefix}.ORF2LCA.txt -o {params.output_prefix}.ORF2LCA.names.txt -t {params.taxonomy_dir}
        singularity exec /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/CAT-5.2/CAT.img CAT add_names -i {params.output_prefix}.ORF2LCA.txt -o {params.output_prefix}.ORF2LCA.officialnames.txt -t {params.taxonomy_dir} --only_official
        singularity exec /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/CAT-5.2/CAT.img CAT summarise -i {params.output_prefix}.contig2classification.officialnames.txt -o {params.output_prefix}.summary.txt -c {input.filtered_scaffolds}
        mv {params.output_prefix}.contig2classification.txt {params.output_prefix}.ORF2LCA.txt {params.output_prefix}.contig2classification.names.txt {params.output_prefix}.contig2classification.officialnames.txt {params.output_prefix}.ORF2LCA.names.txt {params.output_prefix}.ORF2LCA.officialnames.txt {params.output_directory}
        mv {params.output_prefix}.summary.txt {params.output_prefix}.alignment.diamond {params.output_prefix}.log {params.output_prefix}.predicted_proteins.faa {params.output_prefix}.predicted_proteins.gff {params.output_directory}
        """

busco_database_name = basename(normpath(config["busco_database"]))

if config["run_busco"] == True:
    if str(config["busco_version"]) == "3":
        rule busco3:
            input:
                filtered_scaffolds = filtered_scaffolds_filename
            output:
                busco_done = join(assembly_name, "assembly", "busco.done")
            params:
                output_name = "BUSCO3_" + busco_database_name,
                output_directory = "run_" + "BUSCO3_" + busco_database_name,
                dest_directory = join(assembly_name),
                busco_database = config["busco_database"]
            threads: 6
            log: join(assembly_name, "logs", "BUSCO3.log")
            benchmark: join(assembly_name, "benchmarks", "BUSCO3.tsv")
            shell:
                """
                SINGULARITYENV_AUGUSTUS_CONFIG_PATH=/hpc-home/mcgowan/augustus_config singularity exec /hpc-home/mcgowan/software/singularity_testing/busco3.simg run_BUSCO.py -i {input.filtered_scaffolds} -c {threads} -o {params.output_name} -m geno -l {params.busco_database} > {log} 2>&1 && touch {output.busco_done}
                mv {params.output_directory} {params.dest_directory}
                """
    elif str(config["busco_version"]) == "4":
        rule busco4:
            input:
                filtered_scaffolds = filtered_scaffolds_filename
            output:
                busco_done = join(assembly_name, "assembly", "busco.done")
            params:
                output_name = "BUSCO4_" + busco_database_name,
                output_directory = join(assembly_name),
                busco_database = config["busco_database"]
            threads: 6
            log: join(assembly_name, "logs", "BUSCO4.log")
            benchmark: join(assembly_name, "benchmarks", "BUSCO4.tsv")
            shell:
                """
                /ei/projects/e/e5f1ee13-d3bf-4fec-8be8-38c6ad26aac3/data/results/CB-GENANNO-476_DToL_Protists/Software/busco-4.1.2/bin/busco --in {input.filtered_scaffolds} -c {threads} -o {params.output_name} --out_path {params.output_directory} -m geno --offline -l {params.busco_database} > {log} 2>&1 && touch {output.busco_done}
                """
    else:
        print("Error, BUSCO version not supported or recognised:", str(config["busco_version"]))
        sys.exit(-1)

# GENERATE SUMMARY FILE
# For the summary file, show coverage info for each library, and merged libraries if N libraries is > 1
if len(gDNA_samples) > 1:
    gDNA_idxstats_files = (expand(join(assembly_name, "gDNA_alignments", "{sample}.gDNA.idxstats.tsv"), sample = gDNA_samples + ["merged"]))
else:
    gDNA_idxstats_files = (expand(join(assembly_name, "gDNA_alignments", "{sample}.gDNA.idxstats.tsv"), sample = gDNA_samples))

gDNA_depth_file = join(assembly_name, "gDNA_alignments", "depth.tsv")

if len(cDNA_samples) == 0:
    cDNA_idxstats_files = []
    stringtie_assembly_files = []
elif len(cDNA_samples) > 1:
    cDNA_idxstats_files = expand(join(assembly_name, "cDNA_alignments", "{sample}.cDNA.idxstats.tsv"), sample = cDNA_samples + ["merged"])
    stringtie_assembly_files = expand(join(assembly_name, "cDNA_alignments", "{sample}.stringtie.gtf"), sample = cDNA_samples + ["merged"])
else:
    cDNA_idxstats_files = expand(join(assembly_name, "cDNA_alignments", "{sample}.cDNA.idxstats.tsv"), sample = cDNA_samples)
    stringtie_assembly_files = expand(join(assembly_name, "cDNA_alignments", "{sample}.stringtie.gtf"), sample = cDNA_samples)

rule generate_summary:
    input:
        infoseq = join(assembly_name, "assembly", "scaffolds." + str(config["min_scaffold_length"]) + ".tsv"),
        bin_members_parsed = join(assembly_name, "metabat2", "no_coverage", "bin_members.tsv"),
        tiara_classification = join(assembly_name, "tiara", assembly_name + ".tiara.tsv"),
        eukrep_classification = join(assembly_name, "eukrep", assembly_name + ".eukrep.tsv"),
        CAT_classification = join(assembly_name, "CAT", assembly_name + ".CAT.contig2classification.officialnames.txt"),
        blobtools_classification = join(assembly_name, "blobtools", assembly_name + ".blobDB.bestsum.table.txt"), 
        gDNA_depth = gDNA_depth_file,
        gDNA_idxstats = gDNA_idxstats_files,
        cDNA_idxstats = cDNA_idxstats_files,
        stringtie_assemblies = stringtie_assembly_files
    output:
        summary = join(assembly_name, assembly_name + "_summary.tsv")
    script:
        "scripts/generate_summary.py"
