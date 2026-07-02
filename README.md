# nf-RNAseq-pipeline

A Nextflow DSL2 pipeline for bulk RNA-seq preprocessing — from raw paired-end FASTQ files to a gene-level counts matrix ready for differential expression analysis with DESeq2 or edgeR.

---

## Pipeline overview   

```
fastq.gz (paired-end)
    │
    ├─ [0] SAMPLESHEET_CHECK  ──► Auto-detection of R1/R2 pairs
    │
    ├─ [1] FastQC + MultiQC   ──► Pre-trimming QC report
    │
    ├─ [2] TrimGalore          ──► Adapter trimming + quality filtering (Phred ≥ 20)
    │       └─ FastQC          ──► Post-trimming QC report
    │
    ├─ [3] STAR                ──► Splice-aware alignment → sorted BAM
    │
    ├─ [4] SAMtools            ──► BAM indexing + flagstat + idxstats
    │
    ├─ [5] featureCounts       ──► Gene-level quantification → counts matrix
    │
    └─ [6] MultiQC             ──► Aggregated QC report (all steps)

Main outputs:
    results/05_featurecounts/counts_matrix_filtered.txt
    results/06_multiqc/multiqc_report.html
```

Steps [1]–[4] run **in parallel per sample**. Steps [5] and [6] wait for all samples to complete.

---

## Requirements

| Tool | Version | Notes |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.04 | Workflow manager |
| Singularity / Apptainer / Docker | any | Container engine — configurable via profiles |
| Python | ≥ 3.8 | For automatic samplesheet generation |

> All bioinformatics tools (FastQC, TrimGalore, STAR, SAMtools, featureCounts, MultiQC) run inside containers pulled automatically from Quay.io / Seqera — no manual installation required.

> A pre-built STAR genome index and a GTF annotation file are required. See [Building a STAR index](#building-a-star-index) below.

---

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/nmoragas/nf-RNAseq-pipeline.git
cd nf-RNAseq-pipeline
```

### 2. Generate the samplesheet automatically

Point to the directory containing your FASTQ files — R1/R2 pairs are detected automatically when the pipeline runs (no manual samplesheet needed):

```bash
nextflow run main.nf \
    --input /path/to/fastqs/ \
    -profile singularity,slurm
```

If you prefer to inspect the samplesheet before running the full pipeline, you can generate it standalone:

```bash
python bin/generate_samplesheet.py \
    --input /path/to/fastqs/ \
    --output samplesheet.tsv
```

```
sample_id       fastq_r1                        fastq_r2
A549_control    /path/A549_control_R1.fastq.gz  /path/A549_control_R2.fastq.gz
A549_treated    /path/A549_treated_R1.fastq.gz  /path/A549_treated_R2.fastq.gz
```

### 3. Configure reference paths

Edit `nextflow.config` and set your reference files under `params`:

```groovy
params {
    input          = null                                        // set via --input
    outdir         = "${launchDir}/resultats"

    star_genomeDir = "../human_database"                          // pre-built STAR index directory
    counts_gtf     = "../human_database/Homo_sapiens.GRCh38.115.gtf"  // GTF annotation

    fc_strandness  = 0    // 0 = unstranded · 1 = stranded · 2 = reverse-stranded
    min_counts     = 10   // minimum total counts to keep a gene

    trim_cores     = 16
    star_cpus      = 8
    task_cpus      = 8
}
```

### 4. Run the pipeline

```bash
# On a SLURM cluster with Singularity (recommended)
nextflow run main.nf \
    --input /path/to/fastqs/ \
    -profile singularity,slurm

# Resume after interruption (already completed steps are skipped)
nextflow run main.nf \
    --input /path/to/fastqs/ \
    -profile singularity,slurm \
    -resume

# Local execution (for testing with small datasets)
nextflow run main.nf \
    --input /path/to/fastqs/ \
    -profile singularity,standard
```

---

## Repository structure

```
nf-RNAseq-pipeline/
├── main.nf                       ← Main pipeline (Nextflow DSL2)
├── nextflow.config               ← Resources, containers, profiles, reports
├── bin/
│   └── generate_samplesheet.py  ← Auto-generates samplesheet from a FASTQ directory
└── modules/
    ├── sampleSheet_check.nf      ← Samplesheet generation and validation
    ├── fastqc.nf                 ← Read quality control (pre- and post-trim)
    ├── trimgalore.nf             ← Adapter trimming + quality filtering
    ├── star_align.nf             ← Splice-aware genome alignment
    ├── samtools.nf                ← BAM indexing and QC
    ├── featurecounts.nf          ← Gene-level quantification
    └── multiqc.nf                ← Aggregated QC report
```

---

## Output structure

```
resultats/
├── 01_pre_fastqc/
│   └── multiqc_report.html             ← Pre-trimming QC report
├── 02_trimgalore/
│   ├── <sample>_val_1.fq.gz            ← Trimmed R1
│   ├── <sample>_val_2.fq.gz            ← Trimmed R2
│   └── <sample>_trimming_report.txt    ← TrimGalore log
├── 03_star/
│   ├── <sample>.sortedByCoord.out.bam  ← Aligned BAM
│   └── <sample>.Log.final.out          ← STAR alignment summary
├── 04_samtools/
│   ├── <sample>.flagstat               ← Alignment statistics
│   └── <sample>.idxstats               ← Reads per chromosome
├── 05_featurecounts/
│   ├── counts_matrix.txt               ← Full counts matrix (all GTF genes)
│   ├── counts_matrix.txt.summary       ← Read assignment summary
│   └── counts_matrix_filtered.txt      ← Filtered matrix (≥ min_counts) → use this
├── 06_multiqc/
│   └── multiqc_report.html             ← Full pipeline QC report
└── reports/
    ├── execution_report.html           ← Nextflow execution report (resources, status)
    ├── timeline.html                   ← Gantt-style timeline per process
    └── pipeline_trace.txt              ← Per-task resource usage (CPU, memory, time)
```

---

## Key QC metrics to check

After the pipeline completes, open `resultats/06_multiqc/multiqc_report.html` and verify:

| Step | Metric | Acceptable range |
|---|---|---|
| TrimGalore | % reads passing filters | > 95% |
| STAR | % uniquely mapped reads | > 70% |
| featureCounts | % reads assigned to genes | > 60% |
| featureCounts | % multi-mapping reads | < 20% (whole genome) |

---

## Counts matrix format

The file `counts_matrix_filtered.txt` is the direct input for DESeq2 or edgeR in R:

```
Geneid            Length   Sample1   Sample2   Sample3
ENSG00000000003   4535     1523      1891      2034
ENSG00000000005   1687        0         3         1
ENSG00000000419   2356      892      1203       987
```

- **Geneid** — Ensembl gene ID
- **Length** — total exon length in bp (useful for TPM normalisation)
- **Samples** — raw integer read counts

Genes with fewer than `min_counts` total reads across all samples are removed.

---

## Building a STAR index

A STAR genome index is required before running the pipeline. Build it once and reuse across experiments — it is **not** generated by this pipeline.

```bash
# Download genome and annotation from Ensembl
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz

# Build the index (~30 GB RAM · ~1 hour for the human genome)
mkdir -p star_index

STAR \
    --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile Homo_sapiens.GRCh38.115.gtf \
    --sjdbOverhang 149 \
    --runThreadN 8
```

> For a single chromosome (testing only), add `--genomeSAindexNbases 11`

Then point `star_genomeDir` and `counts_gtf` in `nextflow.config` to the generated files.

---

## Execution profiles

This pipeline supports multiple container engines and executors via profiles, combinable with a comma:

```bash
nextflow run main.nf --input /path/to/fastqs/ -profile singularity,slurm
```

| Profile | Effect |
|---|---|
| `standard` | Run locally (`process.executor = 'local'`) |
| `slurm` | Submit each process as a SLURM job (`process.queue = 'odap'`) |
| `singularity` | Use Singularity containers, with bind mounts to `/mnt/typhon/data/references` |
| `apptainer` | Use Apptainer containers (Singularity's successor), same bind mounts |
| `docker` | Docker support (disabled by default — enable in `nextflow.config`) |

---

## Advanced configuration

### Adjust resources per process

Edit the `process` section in `nextflow.config`:

```groovy
withName: 'STAR_ALIGN' {
    cpus   = 16
    memory = '60 GB'
    time   = '8 h'
}
```

### Change the SLURM queue

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue    = 'highmem'   // your queue name
    }
}
```

### Override parameters at runtime

Any parameter in `nextflow.config` can be overridden from the command line without editing the file:

```bash
nextflow run main.nf \
    --input /path/to/fastqs/ \
    --outdir my_results \
    --fc_strandness 1 \
    --min_counts 5 \
    -profile singularity,slurm
```

### Execution reports

Every run automatically generates, under `resultats/reports/`:

- **`execution_report.html`** — per-process CPU, memory and runtime summary
- **`timeline.html`** — Gantt-style visualisation of process execution
- **`pipeline_trace.txt`** — raw per-task resource usage table

Use these to identify bottlenecks and right-size the `cpus`/`memory` settings above.

---

## Tool versions

| Tool | Version | Reference |
|---|---|---|
| FastQC | 0.12.1 | Andrews S. (2010) |
| TrimGalore | 0.6.10 | Krueger F. (2012) |
| STAR | 2.7.11b | Dobin et al. (2013) *Bioinformatics* |
| SAMtools | 1.19 | Li et al. (2009) *Bioinformatics* |
| featureCounts (Subread) | 2.0.6 | Liao et al. (2014) *Bioinformatics* |
| MultiQC | 1.33 | Ewels et al. (2016) *Bioinformatics* |

---

## Citation

If you use this pipeline in your research, please cite the individual tools listed above.

---

## Author

Developed by **Núria Moragas** as part of an RNA-seq bioinformatics training course.
GitHub: [@nmoragas](https://github.com/nmoragas)

Feel free to open an [issue](https://github.com/nmoragas/nf-RNAseq-pipeline/issues) or submit a pull request for suggestions and improvements.
