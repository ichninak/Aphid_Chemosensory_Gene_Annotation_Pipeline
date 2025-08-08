#!/usr/bin/env python3

import os
import glob
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Paramètres par défaut (peuvent être écrasés par config.yaml)
GENOME_DIR = config.get("genome_dir", "/projects/alterevo/00.GROR/genome_data")
SCRIPT_DIR = config.get("script_dir", "/home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/scripts")
DB_CHEMO = config.get("db_chemo", "/home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo")
DBGR = config.get("dbGR", "/home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/dbGR")
DB2GR = config.get("db2GR", "/home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db2GR")
DB2OR = config.get("db2OR", "/home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db2OR")
OUT_BASE = config.get("out_base", "./results")
THREADS = config.get("threads", 50)
THREADS2 = config.get("threads2", 20)

# Flags of execution
RUN_OR = config.get("OR", False)
RUN_GR = config.get("GR", False)

# Guard: GR requires OR
if RUN_GR and not RUN_OR:
    raise ValueError("Config error: GR=true requires OR=true because GR rules depend on OR outputs.")

# Découverte automatique des génomes
genome_files = glob.glob(f"{GENOME_DIR}/*.fa")
genome_files.sort()

# Creation of GAGA identifiers (UN SEUL dictionnaire suffit !)
SAMPLES = {}
for i, genome_path in enumerate(genome_files, 1):
    species = Path(genome_path).stem
    id2 = f"{i:02d}"  # Format 01, 02, 03...
    gaga_id = f"GAGA-00{id2}"  # GAGA-0001, GAGA-0002...
    SAMPLES[gaga_id] = {
        "species": species,
        "genome": genome_path,
        "id": id2
    }

print(f"Discovered {len(SAMPLES)} genomes:")
for gaga_id, info in SAMPLES.items():
    print(f"  {gaga_id}: {info['species']}")

# Helper functions to build family-specific IDs
def gr1_id(sample):
    return f"GAGA-10{SAMPLES[sample]['id']}"

def gr2_id(sample):
    return f"GAGA-20{SAMPLES[sample]['id']}"

def gr3_id(sample):
    return f"GAGA-30{SAMPLES[sample]['id']}"

# Target rules
rule all:
    input:
        # OR workflow results (if enabled)
        expand("{out_base}/{sample}/OR_classification_complete.flag", 
               out_base=OUT_BASE, sample=SAMPLES.keys()) if RUN_OR else [],
        # GR workflow results (if enabled and OR exists)
        expand("{out_base}/{sample}/GR_classification_complete.flag", 
               out_base=OUT_BASE, sample=SAMPLES.keys()) if RUN_GR else []

# Rule OR annotation

# Rule OR annotation - FLAG SIMPLE
rule or_annotation:
    input:
        genome = lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        flag = "{out_base}/{sample}/OR_classification_complete.flag",  # FLAG SIMPLE
        directory = directory("{out_base}/{sample}")
    params:
        gaga_id = lambda wildcards: wildcards.sample,
        species = lambda wildcards: SAMPLES[wildcards.sample]["species"],
        db2or = DB2OR,
        script_dir = SCRIPT_DIR,
        db_chemo = DB_CHEMO,
        out_base = OUT_BASE
    threads: THREADS
    resources:
        mem_mb=20000,
        time_min=1440
    shell:
        """
        set -euo pipefail
        # Run HAPpy ABCENTH for OR
        HAPpy --threads {threads} --annotator ABCENTH --hmm_dir {params.dbor} --genome {input.genome} --output_dir {output.directory}
        
        cd {output.directory}
        
        # Generation of a GFF3 and protein FILE
        echo '----------------------run genometools----------------------'
        gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH.gtf
        gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3
        
        perl {params.script_dir}/gff2fasta_v3.pl {input.genome} ABCENTH_clean.gff3 ABCENTH_clean
        sed 's/X*$//' ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
        mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta
        
        # Validation of OR annotation
        echo '----------------------running Interproscan in the protein set----------------------'
        interproscan.sh -i ABCENTH_clean.pep.fasta -t p -goterms -iprlookup -cpu {threads}
        
        echo '----------------------run BLAST with ORs database----------------------'
        blastp -query ABCENTH_clean.pep.fasta -db {params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads {threads}
        blastp -query ABCENTH_clean.pep.fasta -db {params.db_chemo}/OR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads {threads} -max_target_seqs 5
        blastp -query ABCENTH_clean.pep.fasta -db {params.db_chemo}/GR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads {threads} -max_target_seqs 5
        
        # Classification
        perl {params.script_dir}/run_OR_classification.pl ABCENTH_clean.gff3 {params.gaga_id} {input.genome}
        
        # FLAG SIMPLE : étape terminée 
        touch {output.flag}
        """

# Rules GR1 annotation (ABCENTH)
rule gr1_annotation:
    input:
        or_flag = "{out_base}/{sample}/OR_classification_complete.flag",  # Attend le flag OR
        genome = lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        flag = "{out_base}/{sample}/GR1_classification_complete.flag",
        directory = lambda wildcards: directory(f"{OUT_BASE}/{gr1_id(wildcards.sample)}")
    params:
        gaga_id = lambda wildcards: gr1_id(wildcards.sample),
        species = lambda wildcards: SAMPLES[wildcards.sample]["species"],
        db2gr = DB2GR,
        dbgr = DBGR,
        script_dir = SCRIPT_DIR,
        out_base = OUT_BASE,
        sample = lambda wildcards: wildcards.sample,
        gr1_dir = lambda wildcards: f"{OUT_BASE}/{gr1_id(wildcards.sample)}"
    threads: THREADS
    resources:
        mem_mb=20000,
        time_min=1440
    shell:
        """
        set -euo pipefail
        # Run HAPpy ABCENTH for GR1 (family-specific directory {params.gr1_dir})
        HAPpy --threads {threads} --annotator ABCENTH --hmm_dir {params.db2gr} --genome {input.genome} --output_dir {params.gr1_dir}
        
        cd {params.gr1_dir}

        # Generation of GFF3 and protein files
        gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH.gtf
        gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3
        
        perl {params.script_dir}/gff2fasta_v3.pl {input.genome} ABCENTH_clean.gff3 ABCENTH_clean
        sed 's/X*$//' ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
        mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta
        
        # Classification
        perl {params.script_dir}/run_GR_classification_abcenth_GR_nonamefilter.pl ABCENTH_clean.gff3 {params.gaga_id} {input.genome}
        
        # FLAG SIMPLE : GR1 terminé ✅
        touch {params.out_base}/{params.sample}/GR1_classification_complete.flag
        """

# Rules GR2 annotation (Genewise) 
rule gr2_annotation:
    input:
        gr1_flag = "{out_base}/{sample}/GR1_classification_complete.flag",  # Dépendance GR1 ajoutée
        genome = lambda wildcards: SAMPLES[wildcards.sample]["genome"]
    output:
        flag = "{out_base}/{sample}/GR2_classification_complete.flag",
        directory = lambda wildcards: directory(f"{OUT_BASE}/{gr2_id(wildcards.sample)}")
    params:
        gaga_id = lambda wildcards: gr2_id(wildcards.sample),
        species = lambda wildcards: SAMPLES[wildcards.sample]["species"],
        db2gr = DB2GR,
        dbgr = DBGR,
        script_dir = SCRIPT_DIR,
        out_base = OUT_BASE,
        sample = lambda wildcards: wildcards.sample,
        gr2_dir = lambda wildcards: f"{OUT_BASE}/{gr2_id(wildcards.sample)}"
    threads: THREADS2
    resources:
        mem_mb=10000,
        time_min=1440
    shell:
        """
        set -euo pipefail
        # Run HAPpy Genewise for GR2 (family-specific directory {params.gr2_dir})
        HAPpy --threads {threads} --annotator Genewise --hmm_dir {params.dbgr} --genome {input.genome} --output_dir {params.gr2_dir}
        
        cd {params.gr2_dir}
        
        # Process Genewise results
        perl {params.script_dir}/get_genewise_gtf_corrected_extrafilterfirstexon.pl Genewise.out {input.genome} > Genewise_corrected.gtf
        gt gtf_to_gff3 -o Genewise_corrected.gff3 Genewise_corrected.gtf
        gt gff3 -sort -tidy -retainids -o Genewise_corrected_clean.gff3 Genewise_corrected.gff3
        
        perl {params.script_dir}/get_genewise_cds_to_exon.pl Genewise_corrected_clean.gff3 > Genewise_corrected_clean_cds_to_exon.gff3
        perl {params.script_dir}/gff2fasta_v3.pl {input.genome} Genewise_corrected_clean_cds_to_exon.gff3 Genewise_corrected_clean_cds_to_exon
        
        sed 's/X*$//' Genewise_corrected_clean_cds_to_exon.pep.fasta > Genewise_corrected_clean_cds_to_exon.pep.fasta.tmp
        mv Genewise_corrected_clean_cds_to_exon.pep.fasta.tmp Genewise_corrected_clean_cds_to_exon.pep.fasta
        
        # Classification
        perl {params.script_dir}/run_GR_classification_happy_GR.pl Genewise_corrected_clean_cds_to_exon.gff3 {params.gaga_id} {input.genome}
        
        # FLAG SIMPLE : GR2 terminé ✅
        touch {params.out_base}/{params.sample}/GR2_classification_complete.flag
        """

# Rule GR3 annotation (Combined)
rule gr3_annotation:
    output:
        flag = "{out_base}/{sample}/GR3_classification_complete.flag",
        directory = lambda wildcards: directory(f"{OUT_BASE}/{gr3_id(wildcards.sample)}")
    params:
        gaga_id = lambda wildcards: gr3_id(wildcards.sample),
        species = lambda wildcards: SAMPLES[wildcards.sample]["species"],
        script_dir = SCRIPT_DIR,
        out_base = OUT_BASE,
        sample = lambda wildcards: wildcards.sample,
        gr1_dir = lambda wildcards: f"{OUT_BASE}/{gr1_id(wildcards.sample)}",
        gr2_dir = lambda wildcards: f"{OUT_BASE}/{gr2_id(wildcards.sample)}",
        gr3_dir = lambda wildcards: f"{OUT_BASE}/{gr3_id(wildcards.sample)}"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=60
    shell:
        """
        set -euo pipefail
        mkdir -p {params.gr3_dir}
        cd {params.gr3_dir}
        
        # Combine GR1 and GR2 results (absolute paths)
        perl {params.script_dir}/get_combined_nr_gff.pl {params.gr1_dir}/ABCENTH_clean.gff3 {params.gr2_dir}/Genewise_corrected_clean_cds_to_exon.gff3 > Combined_nr.gff3
        
        perl {params.script_dir}/gff2fasta_v3.pl {input.genome} Combined_nr.gff3 Combined_nr
        sed 's/X*$//' Combined_nr.pep.fasta > Combined_nr.pep.fasta.tmp
        mv Combined_nr.pep.fasta.tmp Combined_nr.pep.fasta
        
        # Get conserved GRs (family ID {params.gaga_id})
        perl {params.script_dir}/get_conserved_GRs.pl Combined_nr.gff3 {params.gaga_id} {input.genome}
        
        touch {output.flag}
        """

# Rules finale pour GR (combine tous les résultats GR)
rule gr_classification:
    input:
        gr1_flag = "{out_base}/{sample}/GR1_classification_complete.flag",
        gr2_flag = "{out_base}/{sample}/GR2_classification_complete.flag",
        gr3_flag = "{out_base}/{sample}/GR3_classification_complete.flag"
    output:
        flag = "{out_base}/{sample}/GR_classification_complete.flag"
    shell:
        """
        set -euo pipefail
        touch {output.flag}
        """
