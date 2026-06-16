/*
============================================================
  Mòdul: FEATURECOUNTS
  Input:  llista de tots els BAMs (collect)
  Output: matriu de counts + summary
============================================================
*/

process FEATURECOUNTS {

    publishDir "${params.outdir}/05_featurecounts", mode: 'copy'

    input:
    path bams      // tots els BAMs junts (col·lectats)
    path gtf

    output:
    path "counts_matrix.txt",         emit: counts
    path "counts_matrix.txt.summary", emit: summary
    path "counts_matrix_filtered.txt",emit: counts_filtered

    script:
    """
    featureCounts \\
        -T ${params.task_cpus} \\
        -p \\
        --countReadPairs \\
        -s ${params.fc_strandness} \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o counts_matrix.txt \\
        ${bams}

    # Neteja la matriu:
    # elimina comentaris, filtra gens amb pocs counts
    # i guarda GeneID + Length + counts
    grep -v "^#" counts_matrix.txt \\
        | awk -v min=${params.min_counts} \\
            'NR==1 { print; next }
             { total=0; for(i=7;i<=NF;i++) total+=\$i;
               if(total>=min) print }' \\
        | cut -f1,6- \\
        > counts_matrix_filtered.txt

    echo "Gens totals (amb capçalera):"
    wc -l counts_matrix.txt

    echo "Gens filtrats (counts >= ${params.min_counts}):"
    wc -l counts_matrix_filtered.txt
    """
}
