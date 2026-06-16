/*
============================================================
  Mòdul: TRIMGALORE
  Input:  tuple (sample_id, [R1, R2])
  Output: reads trimats + reports FastQC post-trim
============================================================
*/

process TRIMGALORE {

    tag "$sample_id"
    publishDir "${params.outdir}/02_trimgalore/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}*_val_1.fq.gz"),
                          path("${sample_id}*_val_2.fq.gz"), emit: reads
    path "*_trimming_report.txt",                                         emit: logs
    path "*_fastqc.{html,zip}",                                           emit: fastqc

    script:
    """
    trim_galore \\
        --paired \\
        --fastqc \\
        --cores ${params.trim_cores} \\
        --quality 20 \\
        --length 20 \\
        ${r1} ${r2}
    """
}

    
