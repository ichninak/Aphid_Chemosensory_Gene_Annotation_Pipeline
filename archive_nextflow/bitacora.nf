nextflow.enable.dsl = 2

workflow bitacora {
    channel
        .fromPath("${params.gene_families_info}/gene_families.xlsx")

}

process bitacoraProcess {
    tag { gagaID1 }
    cpus params.threads
    input:
        tuple val(id2), val(gagaID1), val(species), path(genome), val(PREM)

    output:
        path "${params.out_base}/OBP/${gagaID1}"

    script:
    """
    path="/path/to/Gene-annotation-pipeline"
    table="/path/to/gene_families.xlsx"
    gene_fam_db="/path/to/Gene_families_db"
    genome_name="Genome-name"
    output_directory="${path}/Data/Genomes/Genome-name"
    proteome="${genome_directory}/Proteome-in-fasta-file"
    gff="${genome_directory}/gff3-file"
    genome="${genome_directory}/Genome-in-fasta-file"
    interpro="${genome_directory}/predicted-domains-from-InterPro-in-TSV-format"
    run_bitacora="${path}/bitacora-master/runBITACORA_command_line.sh"
    num_threads=4

    python "$path"/Scripts/run_analysis.py $path $table $gene_fam_db $proteome $interpro $gff $genome $run_bitacora $genome_name $output_directory $num_threads
    """
}