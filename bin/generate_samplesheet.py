#!/usr/bin/env python3
"""
generate_samplesheet.py

Objective:
Scans a directory for fastq.gz files and auto-generates the pipeline CSV sample sheet. Two columns:
    - sample_ID
    - fastq_path

Usage:
    python bin/generate_samplesheet.py --input /path/to/fastqs/
    python bin/generate_samplesheet.py --input /path/to/fastqs/ --output my_samples.tsv
    python bin/generate_samplesheet.py --input /path/to/fastqs/ --single-end
    
    python bin/generate_samplesheet.py --input /data/raw_exemple/ --output samplesheet.tsv
"""



# 0 import librarys
#-------------------------------------------------------------------------
## totes llibraries std de python.
import argparse # Llegir arguments linea comandaments. -- input, -- output, --single-end
import re # Permet treballar amb expressions regulars per cercar i manipular patrons en cadenes de text. Detecció fw / rv. Neteja nom de fitxers.
import sys # Accés a funcions i variables del propi intèrpret de Python. Sortir amb errors (sys.exit) i escriu a stderr
from pathlib import Path # Cerca els FASTQ i manipula rutes


# 1 Patrons d'expressio regulars.
#-------------------------------------------------------------------------
#R1_PAT = re.compile(r'[_\.]R1[_\.]|[_\.]R1\.(fq|fastq)|_1\.(fq|fastq)', re.IGNORECASE)
#R2_PAT = re.compile(r'[_\.]R2[_\.]|[_\.]R2\.(fq|fastq)|_2\.(fq|fastq)', re.IGNORECASE)

R1_PAT = re.compile(r'[_\.\-]R1[_\.\-]|[_\.\-]R1$|[_\.\-]1[_\.\-]|[_\.\-]1$', re.IGNORECASE)
R2_PAT = re.compile(r'[_\.\-]R2[_\.\-]|[_\.\-]R2$|[_\.\-]2[_\.\-]|[_\.\-]2$', re.IGNORECASE)



# 2 Funcions principals:
#-------------------------------------------------------------------------
def find_fastqs(directory):
    """Cerca recursivament tots els fitxers FASTQ."""
    exts = ('*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq')
    files = []
    for ext in exts:
        files.extend(Path(directory).rglob(ext)) 
        # rglob - cerca recursiva (r) a tos els diretoris, subdirctoris. 
        # files.extend — extend acumula resultats a la llista files
    return sorted(set(files)) # set()  elimina duplicats. sorted() odrena alfabeticament


""" Comprovacio:
resultats = find_fastqs('/mnt/hydra/ubs/shared/projects/microbiome/pipelines/shotgun_pipeline/data/raw_data_example')
for i in resultats:
    print(i)
"""

def get_sample_id(filepath, pattern):
    """Elimina el patró R1/R2 i l'extensió per obtenir el sample ID."""
    name = filepath.name # filepath és un objecte Path, i .name extreu només el nom del fitxer sense la ruta:
    # Elimina el patró R1/R2 i tot el que ve després. Usa el patró definit com a separador i es queda amb la part abans del patró ([0]):
    name = pattern.split(name)[0]
    # Elimina guions baixos o punts finals residuals
    name = name.rstrip('_.')
    return name

"""
ex de: 
R1_PAT.split('1000182184C_1.fq.gz')
['1000182184C', None, 'fq', '.gz'] <- es queda amb [0] que es 1000182184C

#prova:
R1_PAT

# Prova amb els teus fitxers reals
fitxer = resultats[0]
print(fitxer)

sid = get_sample_id(fitxer, R1_PAT)
print(f"\nSample ID final: {sid}")
"""


def build_pairs(files):
    """Aparella fitxers R1 amb el seu R2 corresponent."""
    #Pren la llista de tots els fitxers i retorna un diccionari aparellat [nom_mostra : ruta R1 : Ruta R2]
    
    # 1 separar R1 de R2 (a partir llistat directoris de tots els fastq)
    r1s = [f for f in files if R1_PAT.search(f.name)]
    r2s = [f for f in files if R2_PAT.search(f.name)]

    pairs = {}

    # 2 Itera pels R1 i busca el seu R2
    for r1 in r1s:
        sid = get_sample_id(r1, R1_PAT)

        #  3 next() busca el R2 corresponent
        r2_match = next(
            (r2 for r2 in r2s
             if get_sample_id(r2, R2_PAT) == sid
             and r2.parent == r1.parent),
            None
        )

        # 4 Guarda la parella al diccionari
        pairs[sid] = (str(r1.resolve()), str(r2_match.resolve()) if r2_match else None)

    return pairs


""" prova
# resultats: llistat directori dels fastq
r1s = [f for f in resultats if R1_PAT.search(f.name)]
r2s = [f for f in resultats if R2_PAT.search(f.name)]

print("R1s:", r1s)
print("R2s:", r2s)

# Ara construeix les parelles
parelles = build_pairs(resultats)
for sid, (r1, r2) in parelles.items():
    print(f"\nMostra: {sid}")
    print(f"  R1: {r1}")
    print(f"  R2: {r2}")


"""


# 3. Funció `main()` i arguments
#-------------------------------------------------------------------------
def main():

    # 1 Arguments linea de comandaments.
    ## container:
    parser = argparse.ArgumentParser(
        description="Genera un tsv samplesheet per al pipeline de metagenòmica (paired-end)."
    )
    ## Adjunta especificacions
    parser.add_argument('--input',  required=True,             help="Directori amb els fitxers FASTQ")
    parser.add_argument('--output', default='samplesheet.tsv', help="Fitxer TSV de sortida (default: samplesheet.tsv)")

    ## LLegeix l'escrit al terminal, i ho converteix en un objecte python accesible.
        #Quan s'executa:
            # python3 generate_samplesheet.py --input ../data/raw_example --output samplesheet.tsv
            # parse_args() agafa --input ../data/raw_example --output samplesheet.tsv i ho transforma en:
            # args.input   # → '../data/raw_example'
            # args.output  # → 'samplesheet.tsv'
    args = parser.parse_args()

    # 2 Validacion:
    ## El directori existeix?
    if not Path(args.input).is_dir():
        sys.exit(f"ERROR: '{args.input}' no és un directori vàlid.")

    ## Hi ha fitxers FASTQ?
    all_files = find_fastqs(args.input)
    
    if not all_files:
        sys.exit(f"ERROR: No s'han trobat fitxers FASTQ a '{args.input}'.")

    print(f"Trobats {len(all_files)} fitxer(s) FASTQ a {args.input}", file=sys.stderr)



    # 3 Construir parelles R1/R2
    samples = build_pairs(all_files)

    # 4 Construir les files del TSV
    rows = []
    for sid, (r1, r2) in sorted(samples.items()):
        if r2 is None:
            print(f"  AVÍS: No s'ha trobat R2 per '{sid}' — s'omet.", file=sys.stderr)
            continue
        rows.append(f"{sid}\t{r1}\t{r2}\n")

    # 5 Excritura del fitxer:
    with open(args.output, 'w') as fh:
        ## escriptura capçelera
        fh.write("sample_id\tfastq_r1\tfastq_r2\n")
        ## escriu les files 
        fh.writelines(rows)

    print(f"\nSamplesheet escrit → {args.output}  ({len(rows)} mostra/es)")
    # print tots els ID de les mostres:
    #for r in rows:
     #   print(f"  {r.split(chr(9))[0]}")


if __name__ == '__main__':
    main()