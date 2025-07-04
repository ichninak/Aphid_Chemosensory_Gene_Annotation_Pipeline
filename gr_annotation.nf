nextflow.enable.dsl = 2

workflow grWorkflow {
    Channel
        .fromPath("${params.genome_dir}/*.fa")
        .map { it.getName() }
        .sort()
        .map { name -> tuple(name.replaceFirst(/\.fa$/,''), file("${params.genome_dir}/${name}")) }
        .index()
        .map { idx, pair -> tuple(String.format('%02d', idx+1), pair[0], pair[1]) }
        .set { samples }

    samples
        | processHappyGR
        | processAbcenthGR
        | processCombineGR
}

process processHappyGR {
    tag { id }
    cpus params.threads
    input: tuple val(id), val(species), path(genome)
    output: path "${params.out_base}/${id}/01_GR"
    script:
    """
    # ... HAPpy-only GR steps ...
    """
}

process processAbcenthGR {
    tag { id }
    cpus params.threads
    input: tuple val(id), val(species), path(genome)
    output: path "${params.out_base}/${id}/02_GR"
    script:
    """
    # ... ABCENTH-based GR steps ...
    """
}

process processCombineGR {
    tag { id }
    cpus 4
    input: path from processHappyGR, path from processAbcenthGR
    output: path "${params.out_base}/${id}/03_GR"
    script:
    """
    # ... Combine GR results ...
    """
}