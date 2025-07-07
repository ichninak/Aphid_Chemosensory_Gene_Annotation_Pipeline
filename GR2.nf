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
            def gagaID = "GAGA-20${id2}"
            tuple(id4, gagaID, species, genome)
        }
        .set { samples }

    samples
        | processHappyGR
        | processAbcenthGR
        | processCombineGR
}