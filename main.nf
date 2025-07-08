#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { gr1Workflow } from './GR1.nf'
include { gr2Workflow } from './GR2.nf'
include { gr3Workflow } from './GR3.nf'
include { orWorkflow } from './or_annotation.nf'

workflow {
    // run only the steps the user requested
    if( params.GR1 )   gr1Workflow()
    if( params.GR2 )   gr2Workflow()
    if( params.GR3 )   gr3Workflow()
    if( params.OR  )   orWorkflow()

    // if none of option set, send a message to explain how it work
    if ( !params.GR1 && !params.GR2 && !params.GR3 && !params.OR ) {
        log.error 'usage: nextflow run main.nf --option true or false [option: GR1, GR2, GR3 or OR]'
        exit 1
    }
}
