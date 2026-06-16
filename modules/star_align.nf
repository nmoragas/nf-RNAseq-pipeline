/*
============================================================
  Mòdul: STAR_ALIGN
  Input:  tuple (sample_id, R1_trimat, R2_trimat)
  Output: BAM ordenat + log STAR
============================================================
*/

process STAR_ALIGN {

    tag "$sample_id"
    publishDir "${params.outdir}/03_star/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path star_index 
    
    output:
    tuple val(sample_id), path("*.sortedByCoord.out.bam"), emit: bam
    path "*.Log.final.out",                                emit: log
    path "*.SJ.out.tab",                                   emit: sj

    script:
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${star_index} \\
        --readFilesIn ${r1} ${r2} \\
        --readFilesCommand zcat \\
        --runThreadN ${params.star_cpus} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --outFileNamePrefix ${sample_id}. \\
        --outTmpDir tmp_${sample_id}
    """
}

