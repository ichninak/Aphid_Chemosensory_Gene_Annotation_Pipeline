nextflow.enable.dsl = 2

workflow gr1Workflow {
    main:
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
                def PREM = "GAGA-00${id2}"
                tuple(id2, gagaID1, species, genome, PREM)
            }
            .set { samples }

        gr1_out = samples | processAbcenthGR
        
    emit:
        gr1_out
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
    sed s/X*\$// ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
    mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta

    # run interproscan in the protein set

    interproscan.sh -i ABCENTH_clean.pep.fasta -t p -goterms -iprlookup -cpu 40

    # run blast with ORs to obtain names

    blastp -query ABCENTH_clean.pep.fasta -db ${params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
    blastp -query ABCENTH_clean.pep.fasta -db ${params.db_chemo}/GR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5

    cat ${params.db_chemo}/OR_db.fasta ${params.out_base}/${PREM}/${PREM}_ABCENTH_clean_OR_renamed_all.pep.fasta > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query ABCENTH_clean.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query ABCENTH_clean.pep.fasta -db ${params.db_chemo}/GR_dmel_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRdmelblast.txt -num_threads 4 -max_target_seqs 5

    # run script to rename the gff3 and generate the protein file and summary table

    perl ${params.script_dir}/run_GR_classification_abcenth_GR_nonamefilter.pl ABCENTH_clean.gff3 ${PREM} ${genome}
    """
}