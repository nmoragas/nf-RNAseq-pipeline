/*
 * STEP 1/5 — MultiQC
 * Aggregates FastQC reports into a single HTML report.
 * Reusable module: called twice (pre and post trim).
 */
process MULTIQC {

    tag "${step}"

    publishDir "${params.outdir}/${step}/multiqc", mode: 'copy'

    input:
    tuple val(step), path(zips)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data/",       emit: data

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --force
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    """
}
