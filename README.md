# S. suis Five-Antigen Prevalence Pipeline  
**Version v1.0.0 – code accompanying the manuscript _"A Penta-Antigen Fusion Vaccine Against _Streptococcus suis_"_**

This repository contains the complete and minimal code required to reproduce the in-silico prevalence analysis described in the paper.  The workflow searches 5 candidate antigens across high-quality _S. suis_ genomes and summarises their distribution.

## 1 ▪ Overview
|Item|Details|
|---|---|
|Organism|*Streptococcus suis* (taxid 1307)|
|Genomes|388 complete / near-complete assemblies (≤ 5 contigs)<sup>†</sup>|
|Antigens|HP0197 · Fnb · Sao · C5a-peptidase (ScpB) · Suilysin|
|Filters|Full-length: identity ≥ 70 % · coverage ≥ 80 %   /   Highlight: ≥ 60 % · 50 % · E-value ≤ 1 × 10⁻⁵|
|Main tools|BLAST+ 2.15 · Python 3.10 (pandas 2.2.0 · Biopython 1.83) |
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
|`highlight_prevalence_stats.tsv`|Supplementary Data – highlight-domain prevalence (80 % coverage)|

## 4 ▪ Reproducing the Manuscript Results
1. Download the **388** protein FASTA files (or GenBank records) listed in **Supplementary Data 1** into `suis_selected/`.
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

## 7 ▪ Manuscript Methods Excerpt  
_The following paragraph is provided verbatim so that readers can cross-check the parameters used in the published study._

> **Prevalence Analysis of Selected *Streptococcus suis* Antigens**  
> To quantify the prevalence of five candidate antigens (HP0197, Fnbp, Sao, C5a-peptidase, and Suilysin) within a diverse genomic dataset of *S. suis*, we developed a dedicated bioinformatic pipeline, **SSUIS-SADE v1.0.0** (Streptococcus suis Antigen Distribution Explorer).  
>   
> **Dataset curation and assembly.** A comprehensive set of *S. suis* assemblies (taxonomy ID 1307) was downloaded from RefSeq on 1 March 2025 (n = 402). After removing duplicate BioSample IDs and assemblies with a total contig length < 500 kb, 388 genomes (183 serotype-annotated, 205 untyped) remained. The genomic FASTA (*.fna*) records were used as published without further modification.  
>   
> **Sequence database and antigen queries.** The 388 genomic FASTA sequences were concatenated into `all_suis_genomes.fna` and indexed with `makeblastdb` (BLAST+ v2.15.0) [1, 2]. Reference protein sequences for the five full-length antigens were compiled into `query_antigens.fasta`. For a parallel sensitivity analysis, conserved-domain (“highlight”) fragments were compiled into `query_antigens_highlight.fasta`.  
>   
> **Homology search and prevalence metrics.** Protein queries were aligned against the six-frame-translated genome database using `tblastn` (*E*-value ≤ 1 × 10⁻⁵, multi-threaded). Tabular output (format 6) was processed with `parse_prevalence.py` leveraging pandas [3, 4] and Biopython [5]. A genome was considered positive when at least one hit met predefined criteria: full-length ≥ 70 % identity **and** ≥ 80 % coverage; highlight ≥ 60 % identity **and** ≥ 50 % coverage, consistent with previous streptococcal studies [6]. Antigen prevalence was computed as positive genomes / 388 × 100.  
>   
> **Computational reproducibility.** All custom scripts, the Conda environment specification (`environment.yml`), and a Dockerfile (Ubuntu 22.04, BLAST+, Python) are publicly available in the GitHub repository (release v1.0.0). The complete analysis can be executed in one command:  
> ```bash
> docker run --rm -v $(pwd):/work suis-prevalence:1.0 \
>     bash run_suis_prevalence.sh
> ```  
>   
> ### References  
> 1. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. *J Mol Biol.* 1990;215:403-410.  
> 2. Camacho C, Coulouris G, Avagyan V, *et al.* BLAST+: architecture and applications. *BMC Bioinformatics.* 2009;10:421.  
> 3. The pandas development team. *pandas-dev/pandas: Pandas (software).* Zenodo; 2020.  
> 4. McKinney W. Data structures for statistical computing in Python. In: *Proceedings of the 9th Python in Science Conference.* 2010:56-61.  
> 5. Cock PJA, Antao T, Chang JT, *et al.* Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25:1422-1423.  
> 6. Estrada AA, Gottschalk M, Gebhart C, *et al.* Comparative analysis of *Streptococcus suis* genomes identifies novel candidate virulence-associated genes in North American isolates. *Vet Res.* 2022;53:23. https://doi.org/10.1186/s13567-022-01039-8

---
Questions? Open an issue or contact <dlwndghk2056@gmail.com>. 