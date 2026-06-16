# main
// ─── LOG D'INICI ────────────────────────────────────────

/*
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
*/

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
