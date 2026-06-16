/*
============================================================
  Mòdul: SAMTOOLS_QC
  Input:  tuple (sample_id, bam)
  Output: BAM indexat + flagstat + idxstats
============================================================
*/

process SAMTOOLS_QC {

    tag "$sample_id"
    publishDir "${params.outdir}/04_samtools/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("*.bai"), emit: bam_indexed
    path "*.flagstat",                               emit: flagstat
    path "*.idxstats",                               emit: idxstats

    script:
    """
    samtools index -@ ${params.task_cpus} ${bam}

    samtools flagstat \\
        -@ ${params.task_cpus} \\
        ${bam} > ${sample_id}.flagstat

    samtools idxstats ${bam} > ${sample_id}.idxstats
    """
}