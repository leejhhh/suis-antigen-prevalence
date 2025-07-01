# S. suis 5-항원 분포 분석 파이프라인 (논문 부록)

> English version available: [README_EN.md](./README_EN.md)

이 디렉터리는 논문 *「SSUIS Penta-Vaccine MS」* 의 **Methods – In silico analysis** 및 Supplementary Data 재현을 위한 스크립트를 모아 둔 것입니다.

## 📂 핵심 파일
|파일|설명|
|---|---|
|`run_suis_prevalence.sh`|BLAST DB 생성 → tBLASTn 실행 → 결과 파싱까지 자동화한 셸 스크립트|
|`complete_analysis_pipeline.py`|동일 과정을 파이썬 단일 스크립트로 구현(Windows·Linux 호환)|
|`parse_prevalence.py`|BLAST fmt6 결과를 항원별로 필터링·집계|
|`query_antigens.fasta`|5개 항원의 full-length 아미노산 서열|
|`query_antigens_highlight.fasta`|보존 도메인(highlight) 서열|
|`analyze_highlight_sequences.py`|highlight 서열 전용 prevalence 분석|
|`Dockerfile`·`environment.yml`|재현 가능한 실행 환경(Ubuntu 22.04 + Miniconda)|
|`LICENSE`·`CITATION.cff`|MIT 라이선스·인용 메타데이터|

## 🔧 빠른 실행
```bash
# Conda 환경(선택)
conda env create -f environment.yml
conda activate suis_env

# 전체 파이프라인 실행
bash run_suis_prevalence.sh
```
실행 결과는 `suis_prevalence_analysis/` 폴더에 저장됩니다.

## 🗂️ 대용량 데이터
실제 88개 genome FASTA, BLAST 데이터베이스 등(>100 MB)은 GitHub에 포함되지 않았습니다. 

## 📜 인용 방법
이 파이프라인을 사용하실 경우 논문과 Zenodo 릴리스를 함께 인용해 주세요. 자세한 BibTeX은 [README_EN.md](./README_EN.md)를 참고하십시오.

<<<<<<< HEAD
문의: dlwndghk2056@gmail.com
=======
# Expected: prevalence 100 % with 3 genomes listed
```

Pytest will run the same check automatically:

```bash
pytest -q
```

문의: dlwndghk2056@gmail.com
>>>>>>> 52ac12ad80877fc111f3cdf9174fdaf792666ba1
