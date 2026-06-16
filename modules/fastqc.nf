#!/usr/bin/env nextflow

/*
 * STEP 1 / 5 — FastQC
 * Reusable module: called twice (pre-trim and post-trim).
 * The `stage` parameter ("pre_trim" | "post_trim") differentiates output dirs.
 */

process FASTQC {
    // nom logs:
    tag "${step}"
    

    // # Es pot utilitzar tant containers generats: 
    // container "${params.sif_dir}/fastqc.sif" 
    // # com de repositori: 
    //container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"

    //publishDir "${params.outdir}/${step}/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/${step}", mode: 'copy'

    input:
    tuple val(step), val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip

    script:
"""
fastqc \\
    --threads ${task.cpus} \\
    --outdir . \\
    -- ${r1} ${r2}
"""
}