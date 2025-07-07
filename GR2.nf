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
            def gagaID2 = "GAGA-20${id2}"
            def PREM = "GAGA-00${id2}"
            tuple(id4, gagaID2, species, genome, PREM)
        }
        .set { samples }

    samples
        | processHappyGR
}

process processHappyGR {
    tag { gagaID2 }
    cpus params.threads
    input:
        tuple val(id4), val(gagaID2), val(species), path(genome), val(PREM)
    output:
        path "${params.out_base}/${gagaID2}"
    script:
    """
    
    # HAPpy genewise GR

    HAPpy --threads $cpus --genome ${genome} --output_dir ${params.out_base}/${gagaID2} --protein_seqs ${params.db_chemo}/All1_GR.fasta

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