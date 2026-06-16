/*
 * SAMPLESHEET_CHECK
 *
 * Escaneja el directori d'entrada, detecta tots els parells R1/R2
 * independentment de la convenció de noms, i genera un samplesheet
 * TSV validat que el pipeline usarà internament.
 *
 *
 *   Separadors acceptats: _ (guió baix)  . (punt)  - (guió)
 *   Extensions acceptades: .fastq.gz  .fq.gz  .fastq  .fq
 *
 *   R1/R2 al mig del nom (separador abans i després):
 *     sample_R1_001.fastq.gz
 *     sample_S1_L001_R1_001.fastq.gz   (Illumina complet)
 *     sample.R1.001.fastq.gz
 *     sample-R1-001.fastq.gz
 *
 *   R1/R2 al final del nom (separador abans):
 *     sample_R1.fastq.gz
 *     sample_r1.fastq.gz               (insensible a majúscules)
 *     sample.R1.fastq.gz
 *     sample-R1.fastq.gz
 *     sample_R1.fq.gz
 *     sample_R1.fastq
 *
 *   1/2 al mig del nom (separador abans i després):
 *     sample_1_001.fastq.gz
 *     sample.1.extra.fastq.gz
 *     sample-1-x.fastq.gz
 *
 *   1/2 al final del nom (separador abans):
 *     sample_1.fastq.gz
 *     sample.1.fastq.gz
 *     sample-1.fastq.gz
 *     sample_1.fq
 *
 *   NO suportat (ambigu, podria donar falsos positius):
 *     sampleR1.fastq.gz                (sense separador)
 *     sample1.fastq.gz                 (sense separador)
 */

process SAMPLESHEET_CHECK {
    // nom logs:
    tag "00_samplesheet"
    
    
    // Container lleuger amb Python — no necessita eines de bioinformàtica
    // container "docker://python:3.11-slim"
  
    // Publicar el samplesheet generat:
    publishDir "${params.outdir}/00_sample_info", mode: 'copy'

    input:
    val input_dir   // directori amb els fastq.gz

    output:
    path "samplesheet.tsv", emit: samplesheet

    script:
    """
    python3 ${projectDir}/bin/generate_samplesheet.py \\
        --input ${input_dir} \\
        --output samplesheet.tsv
    """
}
