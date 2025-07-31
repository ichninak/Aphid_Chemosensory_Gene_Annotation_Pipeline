nextflow.enable.dsl = 2

workflow gr3Workflow {
    take:
        gr1_out
        gr2_out
    
    main:
        Channel
            .fromPath("${params.genome_dir}/*.fa")
            .map { file -> tuple(file.getName().replaceFirst(/\.fa$/, ''), file) }
            .toSortedList { it[0] }  // Trier par nom de fichier
            .flatMap { list ->
                list.withIndex().collect { item, idx ->
                    def species = item[0]
                    def genome = item[1]
                    def id2 = String.format('%02d', idx+1)
                    def gagaID3 = "GAGA-30${id2}"
                    def PREM = "GAGA-00${id2}"
                    def PREM2 = "GAGA-10${id2}"
                    def PREM3 = "GAGA-20${id2}"
                    tuple(id2, gagaID3, species, genome, PREM, PREM2, PREM3)
                }
            }
            .set { samples }

        gr3_out = samples | processCombineGR
        
    emit:
        gr3_out
}

process processCombineGR {
    tag { gagaID3 }
    cpus params.threads2
    input:
        tuple val(id2), val(gagaID3), val(species), path(genome), val(PREM), val(PREM2), val(PREM3)
    output:
        path "${params.out_base}/${gagaID3}"
    script:
    """
    mkdir ${params.out_base}/${gagaID3}
    cd ${params.out_base}/${gagaID3}

    # Get previous ABCENTH and genewise annotation GFF from the outputs of previous workflows

    cp ${params.out_base}/${PREM2}/${PREM}_ABCENTH_clean_GR_renamed_all_nofragment.gff3 ABCENTH_clean_GR_renamed_all_nofragment.gff3
    cp ${params.out_base}/${PREM2}/${PREM}_GRs.txt ABCENTH_GRs.txt
    cp ${params.out_base}/${PREM3}/${PREM}_Genewise_GR_renamed_all.gff3 Genewise_GR_renamed_all.gff3


    echo "######## Get sugar and conserved GR receptors, and combine with genewise (priority for conserved GRs in ABCENTH)"

    perl ${params.script_dir}/get_conserved_GRs.pl ABCENTH_GRs.txt ABCENTH_clean_GR_renamed_all_nofragment.gff3 

    # It will generate: ABCENTH_GR_all_renamed.gff3 and ABCENTH_GR_Dmelconserved_renamed.gff3

    GeMoMa CompareTranscripts p=Genewise_GR_renamed_all.gff3 a=ABCENTH_GR_Dmelconserved_renamed.gff3 outdir=gemoma_outdir_dmel > gemoma.out 2> gemoma.err

    perl ${params.script_dir}/get_combined_nr_gff.pl Genewise_GR_renamed_all.gff3 ABCENTH_GR_Dmelconserved_renamed.gff3 gemoma_outdir_dmel/comparison.tabular Genewise_AbcenthDmel_combined.gff3


    echo "######## Combine both GFF3 genewise and ABCENTH (priority for gene wise)"

    GeMoMa CompareTranscripts p=ABCENTH_GR_all_renamed.gff3 a=Genewise_AbcenthDmel_combined.gff3 outdir=gemoma_outdir_abcenth > gemoma.out 2> gemoma.err

    perl ${params.script_dir}/get_combined_nr_gff.pl ABCENTH_GR_all_renamed.gff3 Genewise_AbcenthDmel_combined.gff3 gemoma_outdir_abcenth/comparison.tabular GenewiseAbcenth.gff3



    echo "######## Generate the GFF3 and protein file"

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} GenewiseAbcenth.gff3 GenewiseAbcenth
    sed 's/X*\$//' GenewiseAbcenth.pep.fasta > GenewiseAbcenth.pep.fasta.tmp
    mv GenewiseAbcenth.pep.fasta.tmp GenewiseAbcenth.pep.fasta

    echo "######## Run Interpro in the protein set"

    interproscan.sh -i GenewiseAbcenth.pep.fasta -t p -goterms -iprlookup -cpu 4     ## COMMENTED BECAUSE IT IS ALREADY RUN!!!!!


    echo "######## Run blast with ORs to obtain names"


    blastp -query GenewiseAbcenth.pep.fasta -db ${params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.ORcoblast.txt -num_threads 4
    blastp -query GenewiseAbcenth.pep.fasta -db ${params.db_chemo}/GR_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.GRblast.txt -num_threads 4 -max_target_seqs 5

    cat ${params.db_chemo}/OR_db.fasta ${params.out_base}/${PREM}/${PREM}_ABCENTH_clean_OR_renamed_all.pep.fasta > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query GenewiseAbcenth.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query GenewiseAbcenth.pep.fasta -db ${params.db_chemo}/GR_dmel_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.GRdmelblast.txt -num_threads 4 -max_target_seqs 5


    echo "######## Run script to rename the gff3 and generate the protein file and summary table"

    perl ${params.script_dir}/run_GR_classification_happy_GR.pl GenewiseAbcenth.gff3 ${PREM} ${genome}


    """
}