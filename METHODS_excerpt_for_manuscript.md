## Prevalence Analysis of Selected *Streptococcus suis* Antigens across a Strain-Wide Genomic Dataset

To quantify the distribution of the five candidate antigens—HP0197, Fnbp, Sao, C5a-peptidase (ScpB) and Suilysin—in the global *S. suis* population, we implemented a dedicated Python pipeline, **SSUIS-SADE v1.0.0** (Streptococcus suis Antigen Distribution Explorer; GitHub, release v1.0.0).

**Dataset curation.** Complete and near-complete *S. suis* assemblies (taxonomy ID 1307) were downloaded from RefSeq on 1 March 2025 (n = 402). After removing duplicate BioSample IDs and assemblies whose total contig length was < 500 kb, 388 genomes (183 serotype-annotated, 205 untyped) were retained. Each genomic FASTA (*.fna*) record was kept as published; no CDS extraction was performed.

**Database construction.** All 388 genomic sequences were concatenated into a single file (`all_suis_genomes.fna`) and indexed with `makeblastdb` (BLAST+ v2.15.0) [1, 2] using `-dbtype nucl -parse_seqids`.

**Antigen queries.** Reference protein sequences for HP0197 (WP_277937340.1), Fnbp (WP_014636551.1), Sao (WP_211840080.1), ScpB (WP_240208248.1) and Suilysin (AIG43067.1) were compiled in `query_antigens.fasta`.

**tBLASTn search.** Protein queries were aligned against the six-frame-translated genome database with `tblastn` [1, 2]. Searches used an *E*-value threshold of 1 × 10⁻⁵ and eight CPU threads (`-outfmt 6 -num_threads 8`).

**Result parsing and prevalence metrics.** Tabular output (format 6) was processed with the `parse_prevalence.py` module (pandas 2.2.0 [3, 4]; Biopython 1.83 [5]). For each antigen, hits were filtered by identity ≥ 70 % and query-coverage ≥ 80 % (full-length) or identity ≥ 60 % and coverage ≥ 50 % (highlight domains). These thresholds were chosen to balance sensitivity and specificity and mirror those used in previous cross-serotype studies of streptococci [6]. The best HSP per CDS (highest bit-score) was retained, and presence/absence was tallied per genome accession. Prevalence (%) was calculated as (number of genomes with ≥ 1 qualified hit) / 388 × 100. Summary tables for full-length analysis and the 80 %–coverage highlight-domain analysis are provided in Supplementary Data 3 (`highlight_prevalence_stats.tsv`).

**Reproducibility.** All scripts, the `environment.yml`, and a `Dockerfile` (Ubuntu 22.04 base image with BLAST+ and Python dependencies pre-installed) are publicly available in the GitHub release noted above. The entire analysis can be reproduced with:

```bash
docker run --rm -v $(pwd):/work suis-prevalence:1.0 \
    bash run_suis_prevalence.sh
```
requiring only the 388 *.fna* files placed in `suis_selected/`.

---
### References
1. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. *J Mol Biol.* 1990;215(3):403-410. https://doi.org/10.1016/S0022-2836(05)80360-2
2. Camacho C, Coulouris G, Avagyan V, *et al.* BLAST+: architecture and applications. *BMC Bioinformatics.* 2009;10:421. https://doi.org/10.1186/1471-2105-10-421
3. The pandas development team. *pandas-dev/pandas: Pandas (software).* Zenodo; 2020. https://doi.org/10.5281/zenodo.3509134
4. McKinney W. Data structures for statistical computing in Python. In: *Proceedings of the 9th Python in Science Conference.* 2010:56-61. https://doi.org/10.25080/Majora-92bf1922-00a
5. Cock PJA, Antao T, Chang JT, *et al.* Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25(11):1422-1423. https://doi.org/10.1093/bioinformatics/btp163
6. Estrada AA, Gottschalk M, Gebhart C, *et al.* Comparative analysis of *Streptococcus suis* genomes identifies novel candidate virulence-associated genes in North American isolates. *Vet Res.* 2022;53:23. https://doi.org/10.1186/s13567-022-01039-8