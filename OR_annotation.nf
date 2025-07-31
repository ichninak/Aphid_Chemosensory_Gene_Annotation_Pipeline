nextflow.enable.dsl = 2

workflow orWorkflow {
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
                    def gagaID = "GAGA-00${id2}"
                    tuple(id2, gagaID, species, genome)
                }
            }
            .set { samples }

        or_out = samples | processORAnnotation
        
    emit:
        or_out
}

process processORAnnotation {
    tag { gagaID }
    cpus params.threads
    input: 
        tuple val(id2), val(gagaID), val(species), path(genome)
    output: 
        path "${params.out_base}/${gagaID}"
    script:
    """
    # Run HAPpy ABCENTH for OR 

    HAPpy --threads ${task.cpus} --annotator ABCENTH --hmm_dir ${params.dbOR} --genome ${genome} --output_dir ${params.out_base}/${gagaID}

    cd ${params.out_base}/${gagaID}

    # Generation of a GFF3 and protein FILE

    echo '----------------------run genometools----------------------'

    gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH.gtf
    gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} ABCENTH_clean.gff3 ABCENTH_clean
    sed 's/X*\$//' ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
    mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta 

    # Validation of OR annotation

    echo '----------------------running Interproscan in the protein set----------------------'
    
    interproscan.sh -i ABCENTH_clean.pep.fasta -t p -goterms -iprlookup -cpu ${task.cpus}

    echo '----------------------run BLAST with ORs database----------------------'

    blastp -query ${params.out_base}/${gagaID}/ABCENTH_clean.pep.fasta -db ${params.db_chemo}/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads ${task.cpus}
    blastp -query ${params.out_base}/${gagaID}/ABCENTH_clean.pep.fasta -db ${params.db_chemo}/OR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads ${task.cpus} -max_target_seqs 5
    blastp -query ${params.out_base}/${gagaID}/ABCENTH_clean.pep.fasta -db ${params.db_chemo}/GR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads ${task.cpus} -max_target_seqs 5

    # classification 
    
    perl ${params.script_dir}/run_OR_classification.pl ABCENTH_clean.gff3 ${gagaID} ${genome}

    """
}