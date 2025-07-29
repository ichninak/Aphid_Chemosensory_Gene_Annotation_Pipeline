#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { gr1Workflow } from './GR1.nf'
include { gr2Workflow } from './GR2.nf'
include { gr3Workflow } from './GR3.nf'
include { orWorkflow } from './OR_annotation.nf'

// Function to check if OR annotation results exist
def checkORResults() {
    // Check if OR output directories exist with required files
    def baseDir = file("${params.out_base}")
    if (!baseDir.exists()) {
        return false
    }
    
    // Look for GAGA-00XX directories and check for OR results
    def orDirs = baseDir.listFiles().findAll { it.isDirectory() && it.name.matches(/GAGA-00\d+/) }
    
    if (orDirs.isEmpty()) {
        return false
    }
    
    // Check if at least one OR directory has the required output file
    def hasResults = orDirs.any { dir ->
        def requiredFile = file("${dir}/${dir.name}_ABCENTH_clean_OR_renamed_all.pep.fasta")
        if (requiredFile.exists()) {
            log.info "Found OR results in directory: ${dir.name}"
            return true
        }
        return false
    }
    
    return hasResults
}

workflow {
    // run only the steps the user requested
    if( params.GR ) {
        log.info "Starting GR annotation pipeline..."
        
        // Check if OR annotations exist before running GR workflows
        // GR workflows require OR results for BLAST databases
        def orResultsExist = checkORResults()
        
        if (!orResultsExist) {
            log.error """
            ========================================================================================
            ERROR: GR annotation requires OR annotation results to be available first!
            
            OR annotation results not found in output directory: ${params.out_base}
            
            Please run one of the following options first:
            1. Run OR annotation first: nextflow run main.nf --OR true
            2. Run both together: nextflow run main.nf --OR true --GR true
            
            Then you can run GR annotation alone: nextflow run main.nf --GR true
            ========================================================================================
            """.stripIndent()
            exit 1
        }
        
        log.info "OR annotation results found. Proceeding with GR annotation..."
        
        // Execute GR workflows in strict sequence
        log.info "Step 1/3: Running GR1 workflow (ABCENTH annotation)..."
        gr1_result = gr1Workflow()
        
        log.info "Step 2/3: Running GR2 workflow (Genewise annotation)..."
        gr2_result = gr2Workflow(gr1_result)
        
        log.info "Step 3/3: Running GR3 workflow (Combined annotation)..."
        gr3Workflow(gr1_result, gr2_result)
        
        log.info "GR annotation pipeline completed successfully!"
    }
    
    if( params.OR ) {
        log.info "Running OR annotation workflow..."
        or_result = orWorkflow()
        
        // If both OR and GR are requested, run GR after OR completes
        if( params.GR ) {
            log.info "OR completed. Now starting GR annotation pipeline..."
            
            log.info "Step 1/3: Running GR1 workflow (ABCENTH annotation)..."
            gr1_result = gr1Workflow()
            
            log.info "Step 2/3: Running GR2 workflow (Genewise annotation)..."
            gr2_result = gr2Workflow(gr1_result)
            
            log.info "Step 3/3: Running GR3 workflow (Combined annotation)..."
            gr3Workflow(gr1_result, gr2_result)
            
            log.info "Complete OR + GR annotation pipeline finished successfully!"
        }
    }

    // if none of option set, send a message to explain how it works
    if ( !params.GR && !params.OR ) {
        log.error """
        ========================================================================================
        USAGE: nextflow run main.nf [OPTIONS]
        
        OPTIONS:
          --OR true          Run OR annotation workflow
          --GR true          Run GR annotation workflow (requires OR to be done first)
          --OR true --GR true   Run both OR and GR workflows in sequence
          --genome <path>  Path to the genome FASTA file (required for OR and GR workflows)
          --out_base <path>  Base output directory for results (default: ./results)
          --script_dir <path>  Path to the script for OR annotation
          --db_chemo <path>  Path to the chemoreceptor database
          --db2GR <path>  Path to the database for GR annotation
          --dbOR <path>  Path to the database for OR annotation
          --threads <int>  Number of threads to use (default: 50)
          --threads2 <int>  Number of threads to use for second workflow (default: 20)

        EXAMPLES:
          nextflow run main.nf --OR true              # Run OR annotation only
          nextflow run main.nf --GR true              # Run GR annotation (OR must exist)
          nextflow run main.nf --OR true --GR true    # Run both in sequence
        ========================================================================================
        """.stripIndent()
        exit 1
    }
}
