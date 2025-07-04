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
            def id4 = String.format('%04d', idx+1)
            def gagaID = "GAGA-${id4}"
            tuple(id4, gagaID, species, genome)
        }
        .set { samples }

    samples
        | processHappyGR
        | processAbcenthGR
        | processCombineGR
}

process processAbcenthGR {
   tag { gagaID }
    cpus params.threads
    input:
        tuple val(id4), val(gagaID), val(species), path(genome)
    output:
        path "${params.out_base}/${gagaID}"
    script:
    """
    # HAPpy Abcenth GR

    HAPpy --threads $cpus --annotator ABCENTH --hmm_dir ${db2GR} --genome ${genome} --output_dir ${params.out_base}/${gagaID}

    cd ${params.out_base}/${gagaID}

    # Generate a GFF3 and Protein file

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/scripts/get_abcenth_gtf_corrected.pl ABCENTH.gtf $GENOME

    gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH_corrected.gtf
    gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/gff2fasta_v3.pl $GENOME ABCENTH_clean.gff3 ABCENTH_clean
    sed s/X*$// ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
    mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta

    # run interproscan in the protein set

    interproscan.sh -i ABCENTH_clean.pep.fasta -t p -goterms -iprlookup -cpu 40

    # run blast with ORs to obtain names

    blastp -query ABCENTH_clean.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/ORco_sequences.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
    blastp -query ABCENTH_clean.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/All1_GR.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5

    cat /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/All_OR.fasta $ORANNOTATION > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query ABCENTH_clean.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query ABCENTH_clean.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/AGLY_GR.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRaglyblast.txt -num_threads 4 -max_target_seqs 5

    # run script to rename the gff3 and generate the protein file and summary table

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/scripts/run_GR_classification_abcenth_GR_nonamefilter.pl ABCENTH_clean.gff3 $PREM $GENOME
    """
}

process processHappyGR {
    tag { gagaID }
    cpus params.threads
    input:
        tuple val(id4), val(gagaID), val(species), path(genome)
    output:
        path "${params.out_base}/${gagaID}"
    script:
    """
    
    # HAPpy genewise GR

    HAPpy --threads $cpus --genome ${genome} --output_dir ${params.out_base}/${gagaID} --protein_seqs ${params.db_chemo}/All1_GR.fasta

    cd ${params.out_base}/${gagaID}

    # generate a gff3  and protein file 

    cp genewise/all_genewise_predictions.gff .

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/scripts/get_genewise_gtf_corrected_extrafilterfirstexon.pl all_genewise_predictions.gff $GENOME

    gt gtf_to_gff3 -force -o all_genewise_predictions_corrected_unclean.gff3 all_genewise_predictions_corrected.gtf
    gt gff3 -force -sort -tidy -retainids -o all_genewise_predictions_corrected.gff3 all_genewise_predictions_corrected_unclean.gff3

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/gff2fasta_v3.pl $GENOME all_genewise_predictions_corrected.gff3 all_genewise_predictions_corrected
    sed 's/X*$//' all_genewise_predictions_corrected.pep.fasta > all_genewise_predictions_corrected.pep.fasta.tmp
    mv all_genewise_predictions_corrected.pep.fasta.tmp all_genewise_predictions_corrected.pep.fasta

    # convert cds to exon, and get orf to make genes start with M

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/scripts/get_genewise_cds_to_exon.pl all_genewise_predictions_corrected all_genewise_predictions_corrected_cdstoexon.gff3
    gt cds -force -matchdescstart -minorflen 10 -startcodon yes -seqfile $GENOME -o all_genewise_predictions_corrected_cdstoexon_tocds.gff3 all_genewise_predictions_corrected_cdstoexon.gff3
    grep -v 'exon' all_genewise_predictions_corrected_cdstoexon_tocds.gff3 > Genewise.gff3

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/gff2fasta_v3.pl $GENOME Genewise.gff3 Genewise
    sed 's/X*$//' Genewise.pep.fasta > Genewise.pep.fasta.tmp
    mv Genewise.pep.fasta.tmp Genewise.pep.fasta

    # run interproscan in the protein set

    interproscan.sh -i Genewise.pep.fasta -t p -goterms -iprlookup -cpu 40

    # run blast with ORs to obtain names

    blastp -query Genewise.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/ORco_sequences.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORcoblast.txt -num_threads 4
    #blastp -query Genewise.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/OR_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5
    blastp -query Genewise.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/All1_GR.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.GRblast.txt -num_threads 4 -max_target_seqs 5

    cat /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/OR_db.fasta $ORANNOTATION > OR_masAbcenth_db.fasta
    makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
    blastp -query Genewise.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

    blastp -query Genewise.pep.fasta -db /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/db_chemo/AGLY_GR.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.GRaglyblast.txt -num_threads 4 -max_target_seqs 5

    # run script to rename the gff3 and generate the protein and summary table

    perl /home/genouest/inra_umr1349/ichninak/01-GAGA/04_Gene_re-annotation/01-Chemosensory_gene_families/scripts/run_GR_classification_happy_GR.pl Genewise.gff3 $PREM $GENOME

    """
}

process processCombineGR {
    tag { gagaID }
    cpus 4
    input:
      path dir1 from processHappyGR.collect{ it }
      path dir2 from processAbcenthGR.collect{ it }
    output:
      path "${params.out_base}/${gagaID}"
    script:
    """
    mkdir ${params.out_base}/${gagaID}
    cd ${params.out_base}/${gagaID}
    """
}