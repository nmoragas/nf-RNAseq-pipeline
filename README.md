# RNA-seq Preprocessing Pipeline

Pipeline de preprocessament de dades RNA-seq paired-end, des de fitxers FASTQ fins a una matriu de counts llesta per a anàlisi diferencial amb DESeq2.

## Contingut 

```
rnaseq_pipeline/
├── main.nf                  ← Pipeline principal (Nextflow DSL2)
├── nextflow.config          ← Configuració de recursos i contenidors
├── samplesheet.csv          ← Llistat de mostres (generat automàticament)
├── bin/
│   └── generate_samplesheet.py  ← Genera el samplesheet des d'un directori
└── modules/
    ├── fastqc.nf            ← QC inicial dels reads
    ├── trimgalore.nf        ← Trimming d'adaptadors i qualitat
    ├── star.nf              ← Alineament al genoma
    ├── samtools.nf          ← QC del BAM
    ├── featurecounts.nf     ← Quantificació de gens
    └── multiqc.nf           ← Report QC complet
```

---

## Requisits

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Singularity](https://sylabs.io/singularity/) o Docker
- Python >= 3.8 (per generar el samplesheet)
- Índex de STAR ja generat (vegeu secció [Generar índex de STAR](#generar-índex-de-star))

---

## Ús ràpid

### 1. Generar el samplesheet automàticament

Dona el directori on tens els FASTQs i el script detecta automàticament els parells R1/R2:

```bash
python bin/generate_samplesheet.py \
    --input /path/to/fastqs/ \
    --output samplesheet.csv
```

El fitxer generat tindrà aquest format:

```
sample,R1,R2,condition
A549_0_1,/path/A549_0_1_R1.fastq.gz,/path/A549_0_1_R2.fastq.gz,
A549_25_3,/path/A549_25_3_R1.fastq.gz,/path/A549_25_3_R2.fastq.gz,
```

> **Important:** omple la columna `condition` manualment abans de llançar el pipeline.  
> Exemple: `control` per mostres control i `treated` per mostres tractades.

### 2. Configurar els paths de referència

Edita `nextflow.config` i modifica:

```groovy
params {
    star_index    = "/path/to/star_index"      // directori amb l'índex de STAR
    gtf           = "/path/to/annotation.gtf"  // fitxer GTF d'anotació
    outdir        = "results"                  // directori de sortida
    fc_strandness = 0                          // 0=unstranded 1=stranded 2=reverse
}
```

### 3. Llançar el pipeline

```bash
# En local (per proves amb poques mostres)
nextflow run main.nf -profile local

# Al cluster SLURM
nextflow run main.nf -profile slurm

# Reprendre una execució interrompuda (no repeteix passos ja fets)
nextflow run main.nf -profile slurm -resume
```

---

## Flux del pipeline

```
FASTQ (R1 + R2)
    │
    ├─→ FASTQC             QC inicial dels reads crus
    │
    ↓
TRIMGALORE                 Elimina adaptadors i bases de baixa qualitat (Phred < 20)
    │
    ├─→ FASTQC post-trim   QC dels reads nets
    │
    ↓
STAR_ALIGN                 Alineament splice-aware al genoma de referència → BAM
    │
    ↓
SAMTOOLS_QC                Indexació + flagstat + idxstats del BAM
    │
    ↓
FEATURECOUNTS              Compta reads per gen usant el GTF → matriu de counts
    │
    ↓
MULTIQC                    Agrega tots els reports del pipeline en un sol HTML
    │
    ↓
counts_matrix_filtered.txt Matriu final (Gens × Mostres) → input per DESeq2
```

Tots els passos fins a SAMTOOLS_QC s'executen **en paral·lel per mostra**. FEATURECOUNTS i MULTIQC esperen que totes les mostres hagin acabat.

---

## Output

```
results/
├── fastqc/                  QC inicial per mostra
├── trimgalore/              Reads trimats + reports per mostra
├── star/                    BAMs alineats + logs de STAR per mostra
├── samtools/                flagstat + idxstats per mostra
├── featurecounts/
│   ├── counts_matrix.txt           Matriu completa (tots els gens del GTF)
│   ├── counts_matrix.txt.summary   Resum d'assignació de reads
│   └── counts_matrix_filtered.txt  Matriu filtrada (≥ 10 counts) → usar aquesta
└── multiqc/
    └── multiqc_report.html         Report QC complet de tot el pipeline
```

### Mètriques de qualitat a revisar

| Mètrica | Eina | Valor acceptable |
|---|---|---|
| % reads conservats post-trim | TrimGalore | > 95% |
| % uniquely mapped | STAR | > 70% |
| % reads assignats a gens | featureCounts | > 60% |
| % multi-mappers | featureCounts | < 20% (genoma complet) |

---

## Interpretació de la matriu de counts

El fitxer `counts_matrix_filtered.txt` conté:

```
Geneid          Length  Mostra1  Mostra2  ...
ENSG00000000003  4535    1523     1891
ENSG00000000005  1687       0        3
```

- **Geneid**: ID Ensembl del gen
- **Length**: longitud total dels exons en bp (útil per calcular TPM)
- **Mostres**: nombre de reads assignats a cada gen (raw counts enters)

Aquesta matriu és l'input directe per a `DESeq2` o `edgeR` en R.

---

## Generar índex de STAR

L'índex de STAR es genera una sola vegada per genoma i es reutilitza en totes les execucions. No està inclòs al pipeline principal — executa'l per separat:

```bash
# Descarrega el genoma i l'anotació d'Ensembl
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.111.gtf.gz

# Genera l'índex (necessita ~30 GB RAM i ~1h per al genoma humà complet)
mkdir -p star_index

STAR \
    --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile Homo_sapiens.GRCh38.111.gtf \
    --sjdbOverhang 149 \
    --runThreadN 8
```

> Per a un sol cromosoma (proves) afegeix `--genomeSAindexNbases 11`

---

## Configuració avançada

### Canviar recursos per procés

Edita la secció `process` de `nextflow.config`:

```groovy
withName: 'STAR_ALIGN' {
    cpus   = 16       // més CPUs per anar més ràpid
    memory = '60 GB'  // augmenta si el genoma és gran
    time   = '8 h'    // augmenta per a moltes mostres
}
```

### Canviar la cua de SLURM

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue    = 'highmem'   // nom de la teva cua
    }
}
```

### Paràmetres des de la línia de comandes

Pots sobreescriure qualsevol paràmetre sense editar el config:

```bash
nextflow run main.nf \
    -profile slurm \
    --outdir my_results \
    --fc_strandness 1 \
    --min_counts 5
```

---

## Referència de les eines

| Eina | Versió | Referència |
|---|---|---|
| FastQC | 0.12.1 | Andrews S. (2010) |
| TrimGalore | 0.6.10 | Krueger F. (2012) |
| STAR | 2.7.11b | Dobin et al. (2013) Bioinformatics |
| SAMtools | 1.19 | Li et al. (2009) Bioinformatics |
| featureCounts | 2.0.6 | Liao et al. (2014) Bioinformatics |
| MultiQC | 1.25.1 | Ewels et al. (2016) Bioinformatics |
