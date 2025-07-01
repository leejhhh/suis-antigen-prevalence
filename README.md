## S. suis 5-항원 분포 분석 스크립트 (논문 부록)

본 `paper/` 디렉터리는 **「SSUIS Penta-Vaccine MS」** 원고의 *Methods – In silico analysis* 절과 Supplementary Materials 재현을 위해 필요한 최소 파일만을 모아둔 것입니다.

### 1. 파일 설명
|파일|설명|
|----|----|
|`run_suis_prevalence.sh`|BLAST DB 생성, tBLASTn 실행, 결과 TSV 출력까지 수행하는 메인 셸 스크립트|
|`parse_prevalence.py`|BLAST fmt6 결과를 필터링(≥70 % id, ≥80 % cov) 후 항원별 prevalence 계산|
|`complete_analysis_pipeline.py`|위 두 기능을 하나로 묶은 파이썬 버전(선택 사용)|
|`query_antigens.fasta`|5개 항원의 full-length 아미노산 서열|
|`query_antigens_highlight.fasta`|보존 domain(highlight) 서열(선택 분석)|
|`analyze_highlight_sequences.py`|highlight 서열 전용 분석 스크립트|
|`Dockerfile`|Ubuntu 22.04 기반 재현용 컨테이너 정의|
|`README_FULL.md`|프로젝트 전체 설명(원본 README) *참고용*|

### 2. 재현 방법
```bash
# 1) 컨테이너 빌드
cd paper
docker build -t suis_antigen_prevalence .

# 2) 분석 실행 (예: full-length 항원)
docker run --rm -v $(pwd):/app suis_antigen_prevalence bash run_suis_prevalence.sh
```

*로컬 환경(Python 3.10 이상, BLAST+ v2.15+)에서는 `run_suis_prevalence.sh`를 직접 실행해도 됩니다.*

### 3. 논문 첨부 방법
1. 본 `paper/` 폴더를 **ZIP** 으로 압축(`paper.zip`) → *Supplementary Data 1* 으로 업로드
2. 또는 GitHub repository 에 push 후 Release(v1.0.0) → DOI 발급(Zenodo) → 논문 *Data availability* 절에 DOI 기재

### 4. 라이선스 및 인용
- 코드: MIT License (`LICENSE` 참조)
- 인용: `CITATION.cff` 파일 혹은 아래 BibTeX 사용
```bibtex
@software{suis_antigen_prevalence_pipeline,
  author       = {Your Name and Collaborators},
  title        = {S. suis five-antigen prevalence analysis pipeline},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.xxxxxxx}
}
```

### Quickstart (Conda)

```bash
# Clone or unzip the paper/ directory
conda env create -f environment.yml
conda activate suis_env

# Run full pipeline (uses ~4 GB real genomes – ensure data present)
bash run_suis_prevalence.sh
```

### Toy dataset test
A minimal BLAST output file is provided in `sample_data/toy_blast.tsv` so reviewers can verify the parser without downloading any genomes.

```bash
python parse_prevalence.py \
  -i sample_data/toy_blast.tsv \
  -t 3 \
  -q query_antigens.fasta \
  -o sample_data/toy_stats.tsv

# Expected: prevalence 100 % with 3 genomes listed
```

Pytest will run the same check automatically:

```bash
pytest -q
```

문의: dlwndghk2056@gmail.com@gmail.com 
