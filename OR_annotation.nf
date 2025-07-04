nextflow.enable.dsl = 2

workflow orWorkflow {
    Channel
        .fromPath("${params.genome_dir}/*.fa")
        .map { it.getName() }
        .sort()
        .map { name -> tuple(name.replaceFirst(/\.fa$/,''), file("${params.genome_dir}/${name}")) }
        .index()
        .map{ idx, pair ->
            def species = pair[0]
            def genome = pair[1]
            def id4 = String.format('%04d', idx+1)
            def gagaID = "GAGA-${id4}"
            tuple(id4, gagaID, species, genome)
        }
        .set { samples }

    samples | processORAnnotation
}

process processORAnnotation {
    tag { gagaID }
    cpus params.threads
    input: 
        tuple val(id4), val(species), path(genome)
    output: 
        path "${params.out_base}/${gagaID}"
    script:
    """
    # Run HAPpy ABCENTH for OR 

    HAPpy --threads 50 --annotator ABCENTH --hmm_dir ${params.dbOR} --genome ${genome} --output_dir ${params.out_base}/${gagaID}

    cd ${params.out_base}/${gagaID}

    # Generation of a GFF3 and protein FILE

    echo '----------------------run genometools----------------------'

    gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH.gtf
    gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

    perl ${params.script_dir}/gff2fasta_v3.pl ${genome} ABCENTH_clean.gff3 ABCENTH_clean
    sed 's/X*$//' ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
    mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta 

    # Validation of OR annotation

    echo '----------------------running Interproscan in the protein set----------------------'
    
    interproscan.sh -i ABCENTH_clean.pep.fasta -t p -goterms -iprlookup -cpu 40

    echo '----------------------run BLAST with ORs database----------------------'

    blastp -query ${params.out_base}/${gagaID}/ABCENTH_clean.pep.fasta -db ${params.db_chemo}/ORco_sequences.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
    blastp -query ${params.out_base}/${gagaID}/ABCENTH_clean.pep.fasta -db ${params.db_chemo}/OR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 40 -max_target_seqs 5
    blastp -query ${params.out_base}/${gagaID}/ABCENTH_clean.pep.fasta -db ${params.db_chemo}/GR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5

    # classification 
    
    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/run_OR_classification.pl ABCENTH_clean.gff3 ${gagaID} ${genome}

    """
}