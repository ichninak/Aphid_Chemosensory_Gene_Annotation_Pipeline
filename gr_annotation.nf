nextflow.enable.dsl = 2

workflow grWorkflow {
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
            def gagaID1 = "GAGA-10${id2}"
            def gagaID2 = "GAGA-20${id2}"
            def gagaID3 = "GAGA-30${id2}"
            def PREM = "GAGA-00${id2}"
            def PREM2 = "GAGA-10${id2}"
            def PREM3 = "GAGA-20${id2}"
            tuple(id2, gagaID1, gagaID2, gagaID3, species, genome, PREM, PREM2, PREM3)
        }
        .set { samples }

    samples
        | processAbcenthGR
        | processHappyGR
        | processCombineGR
}



process processAbcenthGR {
    tag { gagaID1 }
    cpus params.threads
    input:
        tuple val(id2), val(gagaID1), val(species), path(genome), val(PREM)
    output:
        path "${params.out_base}/${gagaID1}"
    script:
    """
    # HAPpy Abcenth GR

    HAPpy --threads ${task.cpus} --annotator ABCENTH --hmm_dir ${params.db2GR} --genome ${genome} --output_dir ${params.out_base}/${gagaID1}

    cd ${params.out_base}/${gagaID1}

    # Generate a GFF3 and Protein file

    perl ${params.script_dir}/get_abcenth_gtf_corrected.pl ABCENTH.gtf ${genome}

    gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH_corrected.gtf
    gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} ABCENTH_clean.gff3 ABCENTH_clean
    sed s/X*$// ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
    mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta

    # run interproscan in the protein set

    interproscan.sh -i ABCENTH_clean.pep.fasta -t p -goterms -iprlookup -cpu 40

    # run blast with ORs to obtain names

    blastp -query ABCENTH_clean.pep.fasta -db ${params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
    blastp -query ABCENTH_clean.pep.fasta -db ${params.db_chemo}/All1_GR.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5

    cat ${params.db_chemo}/OR_db.fasta ${params.out_base}/${PREM}/${PREM}_ABCENTH_clean_OR_renamed_all.pep.fasta > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query ABCENTH_clean.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query ABCENTH_clean.pep.fasta -db ${params.db_chemo}/AGLY_GR.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRaglyblast.txt -num_threads 4 -max_target_seqs 5

    # run script to rename the gff3 and generate the protein file and summary table

    perl ${params.script_dir}/run_GR_classification_abcenth_GR_nonamefilter.pl ABCENTH_clean.gff3 ${PREM} ${genome}
    """
}


process processHappyGR {
    tag { gagaID2 }
    cpus params.threads
    input:
        tuple val(id2), val(gagaID2), val(species), path(genome), val(PREM)
    output:
        path "${params.out_base}/${gagaID2}"
    script:
    """
    
    # HAPpy genewise GR

    HAPpy --threads ${task.cpus} --genome ${genome} --output_dir ${params.out_base}/${gagaID2} --protein_seqs ${params.db_chemo}/All1_GR.fasta

    cd ${params.out_base}/${gagaID2}

    # generate a gff3  and protein file 

    cp genewise/all_genewise_predictions.gff .

    perl ${params.script_dir}/get_genewise_gtf_corrected_extrafilterfirstexon.pl all_genewise_predictions.gff ${genome}

    gt gtf_to_gff3 -force -o all_genewise_predictions_corrected_unclean.gff3 all_genewise_predictions_corrected.gtf
    gt gff3 -force -sort -tidy -retainids -o all_genewise_predictions_corrected.gff3 all_genewise_predictions_corrected_unclean.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} all_genewise_predictions_corrected.gff3 all_genewise_predictions_corrected
    sed 's/X*$//' all_genewise_predictions_corrected.pep.fasta > all_genewise_predictions_corrected.pep.fasta.tmp
    mv all_genewise_predictions_corrected.pep.fasta.tmp all_genewise_predictions_corrected.pep.fasta

    # convert cds to exon, and get orf to make genes start with M

    perl ${params.script_dir}/get_genewise_cds_to_exon.pl all_genewise_predictions_corrected all_genewise_predictions_corrected_cdstoexon.gff3
    gt cds -force -matchdescstart -minorflen 10 -startcodon yes -seqfile ${genome} -o all_genewise_predictions_corrected_cdstoexon_tocds.gff3 all_genewise_predictions_corrected_cdstoexon.gff3
    grep -v 'exon' all_genewise_predictions_corrected_cdstoexon_tocds.gff3 > Genewise.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} Genewise.gff3 Genewise
    sed 's/X*$//' Genewise.pep.fasta > Genewise.pep.fasta.tmp
    mv Genewise.pep.fasta.tmp Genewise.pep.fasta

    # run interproscan in the protein set

    interproscan.sh -i Genewise.pep.fasta -t p -goterms -iprlookup -cpu 40

    # run blast with ORs to obtain names

    blastp -query Genewise.pep.fasta -db ${params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORcoblast.txt -num_threads 4
    #blastp -query Genewise.pep.fasta -db ${params.db_chemo}/OR_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5
    blastp -query Genewise.pep.fasta -db ${params.db_chemo}/All1_GR.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.GRblast.txt -num_threads 4 -max_target_seqs 5

    cat ${params.db_chemo}/OR_db.fasta ${params.out_base}/${PREM}/${PREM}_ABCENTH_clean_OR_renamed_all.pep.fasta > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query Genewise.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query Genewise.pep.fasta -db ${params.db_chemo}/AGLY_GR.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.GRaglyblast.txt -num_threads 4 -max_target_seqs 5

    # run script to rename the gff3 and generate the protein and summary table

    perl ${params.script_dir}/run_GR_classification_happy_GR.pl Genewise.gff3 ${PREM} ${genome}

    """
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

    # Get previous ABCENTH and genewise annotation GFF


    cp /projects/alterevo/00.GROR/TEST/${PREM2}/${PREM}_ABCENTH_clean_GR_renamed_all_nofragment.gff3 ABCENTH_clean_GR_renamed_all_nofragment.gff3
    cp /projects/alterevo/00.GROR/TEST/${PREM2}/${PREM}_GRs.txt ABCENTH_GRs.txt
    cp /projects/alterevo/00.GROR/TEST/${PREM3}/${PREM}_Genewise_GR_renamed_all.gff3 Genewise_GR_renamed_all.gff3


    echo "######## Get sugar and conserved GR receptors, and combine with genewise (priority for conserved GRs in ABCENTH)"

    perl ${params.script_dir}/get_conserved_GRs.pl ABCENTH_GRs.txt ABCENTH_clean_GR_renamed_all_nofragment.gff3 

    # It will generate: ABCENTH_GR_all_renamed.gff3 and ABCENTH_GR_Aglyconserved_renamed.gff3

    GeMoMa CompareTranscripts p=Genewise_GR_renamed_all.gff3 a=ABCENTH_GR_Aglyconserved_renamed.gff3 outdir=gemoma_outdir_agly > gemoma.out 2> gemoma.err

    perl ${params.script_dir}/get_combined_nr_gff.pl Genewise_GR_renamed_all.gff3 ABCENTH_GR_Aglyconserved_renamed.gff3 gemoma_outdir_agly/comparison.tabular Genewise_AbcenthAgly_combined.gff3


    echo "######## Combine both GFF3 genewise and ABCENTH (priority for gene wise)"

    GeMoMa CompareTranscripts p=ABCENTH_GR_all_renamed.gff3 a=Genewise_AbcenthAgly_combined.gff3 outdir=gemoma_outdir_abcenth > gemoma.out 2> gemoma.err

    perl ${params.script_dir}/get_combined_nr_gff.pl ABCENTH_GR_all_renamed.gff3 Genewise_AbcenthAgly_combined.gff3 gemoma_outdir_abcenth/comparison.tabular GenewiseAbcenth.gff3



    echo "######## Generate the GFF3 and protein file"

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} GenewiseAbcenth.gff3 GenewiseAbcenth
    sed 's/X*$//' GenewiseAbcenth.pep.fasta > GenewiseAbcenth.pep.fasta.tmp
    mv GenewiseAbcenth.pep.fasta.tmp GenewiseAbcenth.pep.fasta

    echo "######## Run Interpro in the protein set"

    interproscan.sh -i GenewiseAbcenth.pep.fasta -t p -goterms -iprlookup -cpu 4     ## COMMENTED BECAUSE IT IS ALREADY RUN!!!!!


    echo "######## Run blast with ORs to obtain names"


    blastp -query GenewiseAbcenth.pep.fasta -db ${params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.ORcoblast.txt -num_threads 4
    blastp -query GenewiseAbcenth.pep.fasta -db ${params.db_chemo}/All1_GR.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.GRblast.txt -num_threads 4 -max_target_seqs 5

    cat ${params.db_chemo}/OR_db.fasta ${params.out_base}/${PREM}/${PREM}_ABCENTH_clean_OR_renamed_all.pep.fasta > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query GenewiseAbcenth.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query GenewiseAbcenth.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/AGLY_GR.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.GRaglyblast.txt -num_threads 4 -max_target_seqs 5


    echo "######## Run script to rename the gff3 and generate the protein file and summary table"

    perl ${params.script_dir}/run_GR_classification_happy_GR.pl GenewiseAbcenth.gff3 ${PREM} ${genome}


    """
}