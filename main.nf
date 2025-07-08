#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { gr1Workflow } from './GR1.nf'
include { gr2Workflow } from './GR2.nf'
include { gr3Workflow } from './GR3.nf'
include { orWorkflow } from './or_annotation.nf'

workflow {
    // run only the steps the user requested
    if( params.GR1 )   happyGR()
    if( params.GR2 )   abcenthGR()
    if( params.GR3 )   combineGR()
    if( params.OR  )   orWorkflow()

    // if none of option set, send a message to explain how it work
    if ( !params.GR1 && !params.GR2 && !params.GR3 && !params.OR ) {
        echo 'usage: XXX'
    }
}
