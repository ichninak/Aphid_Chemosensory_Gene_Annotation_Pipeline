#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { grWorkflow } from './gr_annotation.nf'
include { orWorkflow } from './or_annotation.nf'

workflow {
    // If only OR flag is set, run OR; if only GR flag, run GR; else run both
    if (params.OR && !params.GR) {
        orWorkflow()
    }
    else if (params.GR && !params.OR) {
        grWorkflow()
    }
    else {
        grWorkflow()
        orWorkflow()
    }
}
