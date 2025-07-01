# S. suis Five-Antigen Prevalence Pipeline  
**Version v1.0.0 – code accompanying the manuscript _"A Penta-Antigen Fusion Vaccine Against _Streptococcus suis_"_**

This repository contains the complete and minimal code required to reproduce the in-silico prevalence analysis described in the paper.  The workflow searches 5 candidate antigens across high-quality _S. suis_ genomes and summarises their distribution.

## 1 ▪ Overview
|Item|Details|
|---|---|
|Organism|*Streptococcus suis* (taxid 1307)|
|Genomes|388 complete / near-complete assemblies (≤ 5 contigs)<sup>†</sup>|
|Antigens|HP0197 · Fnb · Sao · C5a-peptidase (ScpB) · Suilysin|
|Filters|identity ≥ 70 % · query-coverage ≥ 80 % · E-value ≤ 1 × 10⁻⁵|
|Main tools|BLAST+ 2.15 · Python 3.10 (pandas · Biopython) |
|Runtime|< 2 min on 8 threads (Intel i7-12700, Ubuntu 22.04)|

† Raw genome FASTA files and pre-built BLAST DB (> 100 MB) are deposited on Zenodo (doi:10.5281/zenodo.XXXXXXX) and are **not** stored in this repository.

## 2 ▪ Quick Start
### 2.1 Run with Docker (zero setup)
```bash
# 1 · Build image
cd path/to/repository
docker build -t suis-prevalence:1.0 .

# 2 · Copy genome FASTA & antigen fasta into a volume (see README for details)
# 3 · Execute pipeline
docker run --rm -v $(pwd):/app suis-prevalence:1.0 bash run_suis_prevalence.sh
```

### 2.2 Run locally
Prerequisites – BLAST+ ≥ v2.15 & Python ≥ 3.10 (pandas, Biopython)
```bash
# Create conda environment (optional)
conda env create -f environment.yml
conda activate suis_env

bash run_suis_prevalence.sh                    # core analysis
python analyze_highlight_sequences.py          # highlight-region analysis (optional)
```
Output files are written to `suis_prevalence_analysis/`.

### 2.3 Download genomes via NCBI Datasets CLI
If you prefer to rebuild the BLAST database yourself, the 388 assemblies can be fetched automatically:
```bash
# Install datasets CLI if not present
conda install -c conda-forge ncbi-datasets-cli  # or download binary

# Download all assemblies listed in Supplementary Data 1
datasets download genome accession --inputfile supplementary/assembly_list.csv --filename suis_genomes.zip
unzip suis_genomes.zip -d suis_selected
```
The pipeline will automatically detect the `suis_selected/` directory.

---
Questions? Open an issue or contact <dlwndghk2056@gmail.com>.

[GitHub release v1.0.0](https://github.com/USER/suis-antigen-prevalence/releases/tag/v1.0.0)

## 3 ▪ Repository Contents
|File / Folder|Purpose|
|-------------|-------|
|`run_suis_prevalence.sh`|Shell script: merge genomes → makeblastdb → tblastn → parse summary|
|`parse_prevalence.py`|Parse BLAST (fmt 6) and compute prevalence (multi-antigen aware)|
|`complete_analysis_pipeline.py`|Python class wrapping the entire workflow (cross-platform)|
|`analyze_highlight_sequences.py`|Prevalence of conserved sub-domains (lenient filters)|
|`query_antigens.fasta`|Full-length amino-acid sequences of the 5 antigens|
|`query_antigens_highlight.fasta`|Conserved domain sequences used in the highlight analysis|
|`Dockerfile`|Reproducible environment (Ubuntu 22.04 + Miniconda + BLAST)|
|`environment.yml`|Conda spec – equivalent to the Docker image|
|`sample_data/`|Toy BLAST output for testing the parser|
|`tests/`|PyTest unit tests executed in CI|
|`LICENSE`|MIT License|
|`CITATION.cff`|Metadata for citation auto-generation|

## 4 ▪ Reproducing the Manuscript Results
1. Download the 88 assembly FASTA files listed in **Supplementary Data 1** into `suis_selected/` (or update the path in the config).
2. Place `query_antigens.fasta` in the repository root (already present).
3. Execute `bash run_suis_prevalence.sh` (or run the Python pipeline).
4. The prevalence table `genomes_with_hit_stats.tsv` and summary report are written to `suis_prevalence_analysis/`.

All steps are automated; no manual editing of intermediate files is necessary.

## 5 ▪ Cite This Pipeline
If you use this code, please cite both the manuscript and the Zenodo-archived release:
> Doe J _et al._ (2025) "Five-Antigen Prevalence Analysis in _Streptococcus suis_." Zenodo. doi:10.5281/zenodo.XXXXXXX.

The `CITATION.cff` file enables automatic citation on GitHub and other services.

## 6 ▪ License
MIT – see `LICENSE` for details. Commercial use is permitted.

---
Questions? Open an issue or contact <dlwndghk2056@gmail.com>. 