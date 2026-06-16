#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*module load apps/nextflow/25.04.6

/*
╔══════════════════════════════════════════════════════════════════════════════╗
║          RNA-seq PREPROCESSING PIPELINE  v1.0                                ║
║                                                                              ║ 
║    FASTQ → counts_matrix                                                     ║ 
║    Paired-end · Nextflow DSL2                                                ║
║                                                                              ║     
║  fastq.gz (paired-end)                                                       ║
║    │                                                                         ║
║    ├─ [0] SAMPLESHEET_CHECK  ──► generate_samplesheet.py                     ║
║    │                              detecció automàtica R1/R2                  ║
║    │                                                                         ║
║    ├─ [1] FastQC  ──► Pre-trim QC                                            ║
║    │                                                                         ║
║    ├─ [2] TrimGalore  ──► Adapter trimming + quality filter (Phred >= 20)    ║
║    │      │                                                                  ║
║    │      └─ [2] FastQC post-trim QC                                         ║
║    │                                                                         ║
║    ├─ [3] STAR  ──► Splice-aware alignment → BAM (sortedByCoord)             ║
║    │                                                                         ║
║    ├─ [4] SAMtools  ──► QC pre qauntification:BAM index + flagstat + idxstats║
║    │                                                                         ║
║    ├─ [5] featureCounts  ──► Gene-level quantification → counts matrix       ║
║    │                                                                         ║
║    └─ [6] MultiQC  ──► Aggregated QC report (all steps)                      ║
║                                                                              ║
║  Output principal:                                                           ║
║    results/05_featurecounts/counts_matrix_filtered.txt                          ║
║    results/06_multiqc_post/multiqc/multiqc_report.html                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
*/


// ─── MODULES IMPORT ────────────────────────────────────────────
include { SAMPLESHEET_CHECK } from './modules/sampleSheet_check.nf'
include { FASTQC as FASTQC_PRE } from './modules/fastqc.nf'
include { MULTIQC as MULTIQC_PRE  } from './modules/multiqc.nf'
include { TRIMGALORE    } from './modules/trimgalore.nf'
include { STAR_ALIGN    } from './modules/star_align.nf'
include { SAMTOOLS_QC   } from './modules/samtools.nf'
include { FEATURECOUNTS } from './modules/featurecounts.nf'
include { FASTQC as FASTQC_POST } from './modules/fastqc.nf'
include { MULTIQC as MULTIQC_POST  } from './modules/multiqc.nf'
/*


*/

// ── Default settings ──────────────────────────────
params.input  = null
params.outdir = "${launchDir}/resultats"

// ── Input validation ──────────────────────────────────────
if (!params.input)      error "ERROR: --input is required (directory with fastq.gz)"



// ─── WORKFLOW  ─────────────────────────────────────────
workflow {

    main:
    // 00. Generar samplesheet automàticament
    ch_input_dir = Channel.fromPath( file(params.input).toAbsolutePath().toString() )
    SAMPLESHEET_CHECK( ch_input_dir )

    
    // Llegir TSV i separar columnes   
    ch_reads = SAMPLESHEET_CHECK.out.samplesheet
    .splitCsv( header: true, sep: '\t' )
    .map { row -> 
        def sample = row.sample_id.replaceAll(/^-/, '')  // elimina - inicial
        [ sample, file(row.fastq_r1), file(row.fastq_r2) ] 
    }

    // 01. FastQC + MULTIQC PRE-alineament   
    FASTQC_PRE( ch_reads.map { sample, r1, r2 -> [ "01_pre_fastqc", sample, r1, r2 ] } )
    MULTIQC_PRE(FASTQC_PRE.out.zip
        .map { sample, zips -> zips }
        .collect()
        .map { zips -> [ "01_pre_fastqc", zips ] }
    )

    
    // 4. Trimming + FastQC post-trim (paral·lel per mostra)
    TRIMGALORE(ch_reads)
    

    // 5. Alineament STAR (paral·lel per mostra)
    //    Combina reads trimats amb l'index de referència
    ch_reads_trimmed = TRIMGALORE.out.reads
    ch_star_index = Channel.fromPath(params.star_genomeDir, type: 'dir', checkIfExists: true)
    STAR_ALIGN(ch_reads_trimmed, ch_star_index.first())

    // 6. QC del BAM (paral·lel per mostra)
    SAMTOOLS_QC(STAR_ALIGN.out.bam)

    // 7. featureCounts — recull TOTS els BAMs i els processa junts
    ch_all_bams = SAMTOOLS_QC.out.bam_indexed
        .map { sample_id, bam, bai -> bam }
        .collect()
    ch_gtf = Channel.fromPath(params.counts_gtf, type: 'dir', checkIfExists: true)
    FEATURECOUNTS(ch_all_bams, ch_gtf.first())



// 8. MultiQC — recull TOTS els reports del pipeline
ch_all_reports = Channel.empty()
    .mix(
        // Extraiem només els zips i els col·lectem
        FASTQC_PRE.out.zip.map { id, zip -> zip }.collect(),
        TRIMGALORE.out.logs.collect(),
        TRIMGALORE.out.fastqc.collect(),
        STAR_ALIGN.out.log.collect(),
        SAMTOOLS_QC.out.flagstat.collect(),
        SAMTOOLS_QC.out.idxstats.collect(),
        FEATURECOUNTS.out.summary.collect()
    )
    .collect() // Ajunta absolutament tots els fitxers en una llista plana
    .map { all_files -> [ "06_multiqc_post", all_files ] } // <--- EL TRUC: Creem la tupla que el teu mòdul espera!

MULTIQC_POST(ch_all_reports)
}




// ─── EVENTS ─────────────────────────────────────────────
workflow.onComplete {
    log.info """
============================================================
  Pipeline completat!
  Estat   : ${workflow.success ? 'OK' : 'ERROR'}
  Output  : ${params.outdir}
  Durada  : ${workflow.duration}
============================================================
""".stripIndent()
}


workflow.onError {
    log.error "Pipeline aturat per error: ${workflow.errorMessage}"
}


// Executar: 
// module load apps/nextflow/25.04.6
// nextflow run ../nf-RNAseq-pipeline/main.nf --input ../fastq_data/ -profile singularity,slurm -resume
// nextflow run ../nf-metagenomics-pipeline/metagenomics.nf --input data/raw_data_example/ -profile apptainer,slurm -resume
