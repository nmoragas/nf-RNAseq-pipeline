#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
============================================================
  RNA-seq Preprocessing Pipeline
  FASTQ → counts_matrix
  Paired-end · Nextflow DSL2
============================================================
*/

// ─── IMPORTS ────────────────────────────────────────────
include { FASTQC        } from './modules/fastqc'
include { TRIMGALORE    } from './modules/trimgalore'
include { STAR_ALIGN    } from './modules/star'
include { SAMTOOLS_QC   } from './modules/samtools'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MULTIQC       } from './modules/multiqc'

// ─── LOG D'INICI ────────────────────────────────────────
log.info """
============================================================
  RNA-seq Pipeline
============================================================
  Samplesheet : ${params.samplesheet}
  STAR index  : ${params.star_index}
  GTF         : ${params.gtf}
  Output dir  : ${params.outdir}
  Strandness  : ${params.fc_strandness}
  Min counts  : ${params.min_counts}
============================================================
""".stripIndent()

// ─── WORKFLOW PRINCIPAL ──────────────────────────────────
workflow {

    // 1. Llegeix el samplesheet i crea el canal de mostres
    //    Format CSV: sample,R1,R2,condition
    ch_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def r1        = file(row.R1)
            def r2        = file(row.R2)
            return [ sample_id, [ r1, r2 ] ]
        }

    // Valida que els fitxers existeixen
    ch_samples.view { sample_id, reads ->
        "Mostra: ${sample_id} | R1: ${reads[0].name} | R2: ${reads[1].name}"
    }

    // 2. Referència
    ch_star_index = Channel.fromPath(params.star_index, type: 'dir')
    ch_gtf        = Channel.fromPath(params.gtf)

    // 3. FastQC pre-trim (paral·lel per mostra)
    FASTQC(ch_samples)

    // 4. Trimming + FastQC post-trim (paral·lel per mostra)
    TRIMGALORE(ch_samples)

    // 5. Alineament STAR (paral·lel per mostra)
    //    Combina reads trimats amb l'index de referència
    ch_reads_trimmed = TRIMGALORE.out.reads
    STAR_ALIGN(ch_reads_trimmed, ch_star_index.first())

    // 6. QC del BAM (paral·lel per mostra)
    SAMTOOLS_QC(STAR_ALIGN.out.bam)

    // 7. featureCounts — recull TOTS els BAMs i els processa junts
    ch_all_bams = SAMTOOLS_QC.out.bam_indexed
        .map { sample_id, bam, bai -> bam }
        .collect()

    FEATURECOUNTS(ch_all_bams, ch_gtf.first())

    // 8. MultiQC — recull TOTS els reports del pipeline
    ch_all_reports = Channel.empty()
        .mix(
            FASTQC.out.zip.map        { id, zip  -> zip  }.collect(),
            TRIMGALORE.out.logs.collect(),
            TRIMGALORE.out.fastqc.collect(),
            STAR_ALIGN.out.log.collect(),
            SAMTOOLS_QC.out.flagstat.collect(),
            SAMTOOLS_QC.out.idxstats.collect(),
            FEATURECOUNTS.out.summary.collect()
        )
        .collect()

    MULTIQC(ch_all_reports)
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
