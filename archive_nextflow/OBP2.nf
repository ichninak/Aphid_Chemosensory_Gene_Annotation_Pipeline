nextflow.enable.dsl = 2

workflow obpWorkflow {
    Channel
        .fromPath("${params.genome_dir}/*.fa")
        .map { it.getName() }
        .sort()
        .map { name -> tuple(name.replaceFirst(/\.fa$/,''), file("${params.genome_dir}/${name}")) }
        .index()
        .map { idx, pair ->
            def species = pair[0]
            def genome = pair[1]
            def id2 = String.format('%02d', idx+1)
            def gagaID1 = "GAGA-00${id2}"
            def PREM = "GAGA-00${id2}"
            tuple(id2, gagaID1, species, genome, PREM)
        }
        .set { samples }

    samples | processObp
}

process processObp {
    tag { gagaID1 }
    cpus params.threads
    input:
        tuple val(id2), val(gagaID1), val(species), path(genome), val(PREM)

    output:
        path "${params.out_base}/OBP/${gagaID1}"

    script:
    """

    mkdir -p ${params.out_base}/OBP/${gagaID1}
    cd ${params.out_base}/OBP/${gagaID1}

    # Get de annotation from BITACORA already run in the re-annotation pipeline

    cp /path/to/re-annotation/OBP_genomic_and_annotated_proteins_trimmed_idseqsclustered.gff3 Bitacora_raw.gff3

    # Using untrimmed proteins
    #cat /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/${gagaID1}/gene_families_pipeline/OBP/Step2_bitacora/OBP/Intermediate_files/OBP_annot_genes.gff3 /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/${gagaID1}}/gene_families_pipeline/OBP/Step2_bitacora/OBP/Intermediate_files/OBP_genomic_genes.gff3 > Bitacora_raw.gff3

    ######## Generate a GFF3 and protein file ########

    sed s/split/separated/g Bitacora_raw.gff3 > Bitacora_raws.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} Bitacora_raws.gff3 Bitacora_raws
    sed s/X*$// Bitacora_raws.pep.fasta > Bitacora_raws.pep.fasta.tmp
    mv Bitacora_raws.pep.fasta.tmp Bitacora_raws.pep.fasta

    # identify and curate chimeric sequences
    
    blastp -subject Bitacora_raws.pep.fasta -query ${params.db_chemo}/OBP_db.fasta -out Bitacora_raws.pep.fasta.iblast.txt -outfmt "6 std qlen slen" -evalue 1e-5

    #perl /home/projects/ku_00039/people/joeviz/programs/bitacora_modftrimlength/Scripts/get_blastp_parsed_newv2.pl Bitacora_raws.pep.fasta.iblast.txt Bitacora_raws.pep.fasta.iblast 1e-5
    #perl /home/projects/ku_00039/people/joeviz/programs/bitacora_modftrimlength/Scripts/get_blast_hmmer_combined.pl Bitacora_raws.pep.fasta.iblastblastp_parsed_list.txt Bitacora_raws.pep.fasta.iblast

    perl ${params.script_dir}/get_fullproteinlist_curatingchimeras_trimsingle.pl Bitacora_raws.pep.fasta Bitacora_raws.pep.fasta.iblast_combinedsearches_list.txt
    perl ${params.script_dir}/get_annot_genes_gff_v2.pl Bitacora_raws.gff3 ${genome} Bitacora_raws.pep.fasta.iblast_combinedsearches_full_fixed_list.txt Bitacora
    
    sed s/separated/split/g Bitacora_annot_genes_trimmed.gff3 > Bitacora.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} Bitacora.gff3 Bitacora
    sed s/X*$// Bitacora.pep.fasta > Bitacora.pep.fasta.tmp
    mv Bitacora.pep.fasta.tmp Bitacora.pep.fasta

    ######## Run interproscan in the protein set ########
    
    interproscan.sh -i Bitacora.pep.fasta -t p -goterms -iprlookup -cpu 4

    ######## Run blast with ORs to obtain names ########

    blastp -query Bitacora.pep.fasta -db ${params.db_chemo}/OBP_db.fasta -outfmt "6 std qlen slen" -out Bitacora.pep.fasta.OBPblast.txt -num_threads 4
    blastp -query Bitacora.pep.fasta -db ${params.db_chemo}/OBP_incorrect_venom_allergen_db.fasta -outfmt "6 std qlen slen" -out Bitacora.pep.fasta.OBPvenomallergenblast.txt -num_threads 4

    ######## Run script to rename the gff3 and generate the protein file and summary table ########

    perl ${params.script_dir}/run_classification_bitacora_OBP.pl Bitacora.gff3 ${gagaID1} ${genome}

    """
}