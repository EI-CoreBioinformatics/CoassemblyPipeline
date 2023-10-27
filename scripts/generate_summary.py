#!/usr/bin/env python

# This script generates a summary file from the co-assembly pipeline
# gDNA and cDNA alignment stats to be reported last as the number of columns depends on number of samples

import os
import sys

output_file = snakemake.output.summary
fo = open(output_file, "w")

contigs = {}
all_contigs = []
gDNA_samples = []
cDNA_samples = []

# Length and GC content, list of all contigs
with open(snakemake.input.infoseq, "r") as f:
    f.readline()

    for line in f:
        line = line.strip().split()

        seq_id = line[0]
        seq_length = line[1]
        seq_gc = line[2]

        all_contigs.append(seq_id)

        contigs[seq_id] = {}
        contigs[seq_id]["length"] = seq_length
        contigs[seq_id]["gc"] = seq_gc

# MetaBat2 bin membership
with open(snakemake.input.bin_members_parsed, "r") as f:
    f.readline()

    for line in f:
        line = line.strip().split("\t")

        seq_id = line[0]
        bin_member = line[1]

        contigs[seq_id]["bin"] = bin_member

# Tiara classification
# Special case, if classified as organellar, also input type (class_snd_stage)
with open(snakemake.input.tiara_classification, "r") as f:
    f.readline()
    
    for line in f:
        line = line.strip().split("\t")

        seq_id = line[0]
        tiara_classification = line[1]
        tiara_second_stage = line[2]

        if tiara_classification == "organelle":
            tiara_classification = tiara_classification + "_" + tiara_second_stage

        contigs[seq_id]["tiara"] = tiara_classification

# EukRep classification
with open(snakemake.input.eukrep_classification, "r") as f:
    f.readline()
    
    for line in f:
        line = line.strip().split("\t")

        seq_id = line[0]
        eukrep_classification = line[1]

        contigs[seq_id]["eukrep"] = eukrep_classification

# CAT classification
with open(snakemake.input.CAT_classification, "r") as f:
    f.readline()

    for line in f:
        line = line.strip().split("\t")

        seq_id = line[0]
        classification = line[1]
        cat_reason = line[2]

        if classification == "no taxid assigned":
            cat_superkingom = "Unclassified"
            cat_phylum = "Unclassified"
            cat_class = "Unclassified"
            cat_order = "Unclassified"
            cat_family = "Unclassified"
            cat_genus = "Unclassified"
            cat_species = "Unclassified"
        else:
            cat_superkingom = line[5]
            if cat_superkingom != "no support":
                cat_superkingom = cat_superkingom.split(":")[0]
            
            cat_phylum = line[6]
            if cat_phylum != "no support":
                cat_phylum = cat_phylum.split(":")[0]
            
            cat_class = line[7]
            if cat_class != "no support":
                cat_class = cat_class.split(":")[0]
            
            cat_order = line[8]
            if cat_order != "no support":
                cat_order = cat_order.split(":")[0]

            cat_family = line[9]
            if cat_family != "no support":
                cat_family = cat_family.split(":")[0]
            
            cat_genus = line[10]
            if cat_genus != "no support":
                cat_genus = cat_genus.split(":")[0]

            cat_species = line[11]
            if cat_species != "no support":
                cat_species = cat_species.split(":")[0]

        contigs[seq_id]["cat_superkindom"] = cat_superkingom
        contigs[seq_id]["cat_phylum"] = cat_phylum
        contigs[seq_id]["cat_class"] = cat_class
        contigs[seq_id]["cat_order"] = cat_order
        contigs[seq_id]["cat_family"] = cat_family
        contigs[seq_id]["cat_genus"] = cat_genus
        contigs[seq_id]["cat_species"] = cat_species
        contigs[seq_id]["cat_reason"] = cat_reason

# Blobtools classification
# First 11 lines are comments, 12th is header
with open(snakemake.input.blobtools_classification, "r") as f:
    for i in range(12):
        f.readline()
    
    for line in f:
        line = line.strip().split("\t")

        seq_id = line[0]
        blobtools_superkingdom = line[5]
        blobtools_phylum = line[8]
        blobtools_order = line[11]
        blobtools_family = line[14]
        blobtools_genus = line[17]
        blobtools_species = line[20]

        contigs[seq_id]["blobtools_superkingdom"] = blobtools_superkingdom
        contigs[seq_id]["blobtools_phylum"] = blobtools_phylum
        contigs[seq_id]["blobtools_order"] = blobtools_order
        contigs[seq_id]["blobtools_family"] = blobtools_family
        contigs[seq_id]["blobtools_genus"] = blobtools_genus
        contigs[seq_id]["blobtools_species"] = blobtools_species

# gDNA coverage depth
gDNA_depth = {} # have a separate dict to store 
gDNA_depth_header = ""
with open(snakemake.input.gDNA_depth, "r") as f:
    gDNA_depth_header = f.readline().strip().split()[2:]

    for line in f:
        line = line.strip().split()

        seq_id = line[0]
        
        line = line[2:]

        gDNA_depth[seq_id] = line

# gDNA idxstats per sample
for idxstats in snakemake.input.gDNA_idxstats:
    sample_name = os.path.basename(idxstats).split(".gDNA.idxstats.tsv")[0]
    gDNA_samples.append(sample_name)

    with open(idxstats, "r") as f:
        for line in f:
            line = line.strip().split("\t")

            seq_id = line[0]
            alignments = line[2]

            if seq_id != "*":
                contigs[seq_id]["gDNA_idxstats_" + sample_name] = alignments

# cDNA idxstats per sample
for idxstats in snakemake.input.cDNA_idxstats:
    sample_name = os.path.basename(idxstats).split(".cDNA.idxstats.tsv")[0]
    cDNA_samples.append(sample_name)

    with open(idxstats, "r") as f:
        for line in f:
            line = line.strip().split("\t")

            seq_id = line[0]
            alignments = line[2]

            if seq_id != "*":
                contigs[seq_id]["cDNA_idxstats_" + sample_name] = alignments

# Stringtie transcripts per sample per contig
for stringtie_assembly in snakemake.input.stringtie_assemblies:
    sample_name = os.path.basename(stringtie_assembly).split(".stringtie.gtf")[0]

    transcripts_per_contig = {}
    with open(stringtie_assembly, "r") as f:
        for line in f:
            if line[0] != "#":
                line = line.strip().split()

                seq_id = line[0]
                feature_type = line[2]

                if feature_type == "transcript":
                    if seq_id in transcripts_per_contig:
                        transcripts_per_contig[seq_id] += 1
                    else:
                        transcripts_per_contig[seq_id] = 1

        for seq_id in all_contigs:
            if seq_id in transcripts_per_contig:
                contigs[seq_id]["stringtie_" + sample_name] = transcripts_per_contig[seq_id]
            else:
                contigs[seq_id]["stringtie_" + sample_name] = 0

header = ["sequence_id", "length", "GC", "bin", "tiara", "eukrep", "cat_superkingdom", "cat_phylum", "cat_class"]
header += ["cat_order", "cat_family", "cat_genus", "cat_species", "cat_reason", "blobtools_superkingdom"]
header += ["blobtools_phylum", "blobtools_order", "blobtools_family", "blobtools_genus", "blobtools_species"]

for sample_name in cDNA_samples:
    header += ["cDNA_idxstats_" + sample_name]
    header += ["stringtie_transcripts_" + sample_name]

for sample_name in gDNA_samples:
    header += ["gDNA_idxstats_" + sample_name]

for h in gDNA_depth_header:
    header += ["depth_" + h]

fo.write("\t".join(header) + "\n")

for seq_id in all_contigs:
    line = [seq_id]

    line.append(contigs[seq_id]["length"])
    line.append(contigs[seq_id]["gc"])
    line.append(contigs[seq_id]["bin"])
    line.append(contigs[seq_id]["tiara"])
    line.append(contigs[seq_id]["eukrep"])
    line.append(contigs[seq_id]["cat_superkindom"])
    line.append(contigs[seq_id]["cat_phylum"])
    line.append(contigs[seq_id]["cat_class"])
    line.append(contigs[seq_id]["cat_order"])
    line.append(contigs[seq_id]["cat_family"])
    line.append(contigs[seq_id]["cat_genus"])
    line.append(contigs[seq_id]["cat_species"])
    line.append(contigs[seq_id]["cat_reason"])
    line.append(contigs[seq_id]["blobtools_superkingdom"])
    line.append(contigs[seq_id]["blobtools_phylum"])
    line.append(contigs[seq_id]["blobtools_order"])
    line.append(contigs[seq_id]["blobtools_family"])
    line.append(contigs[seq_id]["blobtools_genus"])
    line.append(contigs[seq_id]["blobtools_species"])

    for sample_name in cDNA_samples:
        line.append(contigs[seq_id]["cDNA_idxstats_" + sample_name])
        line.append(contigs[seq_id]["stringtie_" + sample_name])

    for sample_name in gDNA_samples:
        line.append(contigs[seq_id]["gDNA_idxstats_" + sample_name])

    line += gDNA_depth[seq_id]

    fo.write("\t".join(map(str, line)) + "\n")
