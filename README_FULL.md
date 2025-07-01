# S. suis 5-Antigen Prevalence Analysis

이 프로젝트는 *Streptococcus suis* 5개 주요 항원의 분포를 분석하는 파이프라인입니다.

## 📊 프로젝트 개요

**목표**: 5개의 주요 S. suis 항원이 Complete genome + 5개 이하 contig를 가진 88개 고품질 genome에서 얼마나 분포하고 있는지 %로 분석

**분석 대상 항원**:
- HP0197 (WP_277937340.1)
- Fnb (WP_014636551.1) 
- SAO (WP_211840080.1)
- C5a (WP_240208248.1)
- Suilysin (AIG43067.1)

**분석 기준**:
- Identity ≥ 70%
- Coverage ≥ 80%
- E-value ≤ 1e-5

## 📁 프로젝트 구조

1.  **데이터 준비 (완료):**
    *   NCBI Datasets에서 S. suis 메타데이터 수집
    *   88개 고품질 조립체 선별 (Complete Genome + contig ≤ 5)
    *   Genome 조립체 (`.fna` 파일) 다운로드 완료
2.  **분석 파이프라인:**
    *   88개 `.fna` 파일을 단일 FASTA로 병합 (`all_suis_genomes.fna`)
    *   BLAST 데이터베이스 생성
    *   `query_antigens.fasta`를 이용한 `tblastn` 검색
    *   `parse_prevalence_by_antigen.py`로 항원별 분포율 계산

## 🚀 빠른 실행 방법

### Docker를 이용한 분석 실행

```bash
# 1. Docker 컨테이너 실행
docker run --rm -it -v "${PWD}:/app" suis_prevalence_analyzer bash

# 2. 컨테이너 내에서 분석 실행  
bash run_antigen_analysis.sh
```

### 개별 스크립트 실행

```bash
# 기본 prevalence 분석
bash run_suis_prevalence.sh

# 항원별 상세 분석
bash run_antigen_analysis.sh
```

## 📂 파일 구조

```
offline_bundle/
├─ Miniconda3-latest-Linux-x86_64.sh  # Miniconda installer for Linux
├─ ncbi-blast-2.15.0+-x64-linux.tar.gz # BLAST+ binaries for Linux
├─ datasets-linux-amd64              # NCBI Datasets CLI executable for Linux
├─ run_suis_prevalence.sh            # Main execution script
├─ parse_prevalence.py               # Python script for parsing BLAST results
├─ query_antigens.fasta              # FASTA file with 5 query antigen sequences
└─ suis_selected/                    # Directory containing 388 .fna genome files
    ├─ GCF_000000001.1.fna
    └─ ... (all 388 files)
```

**Action Required:** Copy your project-specific files into the `C:\Users\LJH\suis_prevalence\offline_bundle` directory:

1.  Copy your `run_suis_prevalence.sh` script.
2.  Copy your `parse_prevalence.py` script.
3.  Copy your `query_antigens.fasta` file.
4.  Copy your entire `suis_selected/` directory (containing the 388 `.fna` files).

## Offline Linux Machine Setup and Execution

Once the `offline_bundle` directory is prepared, transfer it to the target offline Linux machine.

**Option 1: Manual Installation and Execution**

1.  **Navigate to the bundle directory:**
    ```bash
    cd /path/to/offline_bundle
    ```

2.  **Install Miniconda:**
    ```bash
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    # Add conda to PATH (adjust if you installed elsewhere or prefer not to modify .bashrc)
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    # Or for the current session only:
    # export PATH="$HOME/miniconda3/bin:$PATH"
    ```

3.  **Install Required System Packages (if not present):**
    You might need `libgomp1`. Use the system's package manager (e.g., `apt` for Debian/Ubuntu, `yum`/`dnf` for CentOS/Fedora). This step might require temporary internet access or offline package repositories if the base OS doesn't include it.
    ```bash
    # Example for Ubuntu/Debian (requires sudo/root)
    # sudo apt-get update && sudo apt-get install -y libgomp1
    ```

4.  **Install Python Libraries:**
    ```bash
    pip install pandas biopython
    ```
    *(Note: If the offline machine has never had pip installed or lacks build tools, installing wheels directly might be necessary. You could download the `.whl` files for pandas, biopython, and their dependencies on an online machine using `pip download pandas biopython --platform manylinux_2_17_x86_64 --only-binary=:all:` and transfer them.)*

5.  **Install BLAST+:**
    ```bash
    mkdir -p $HOME/tools/blast
    tar xzf ncbi-blast-2.15.0+-x64-linux.tar.gz -C $HOME/tools/blast --strip-components=1
    # Add BLAST to PATH
    echo 'export PATH="$HOME/tools/blast/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    # Or for the current session only:
    # export PATH="$HOME/tools/blast/bin:$PATH"
    ```

6.  **Install NCBI Datasets CLI:**
    ```bash
    mkdir -p $HOME/tools/ncbi-datasets
    cp datasets-linux-amd64 $HOME/tools/ncbi-datasets/datasets
    chmod +x $HOME/tools/ncbi-datasets/datasets
    # Add datasets to PATH
    echo 'export PATH="$HOME/tools/ncbi-datasets:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    # Or for the current session only:
    # export PATH="$HOME/tools/ncbi-datasets:$PATH"
    ```

7.  **Make Scripts Executable:**
    ```bash
    chmod +x run_suis_prevalence.sh parse_prevalence.py
    ```

8.  **Run the Analysis:**
    ```bash
    bash run_suis_prevalence.sh
    ```
    The script will perform the FASTA concatenation, BLAST DB creation, `tblastn` search, and prevalence calculation. Results will be printed to standard output or saved to files as defined in the script.

**Option 2: Docker-Based Execution**

This method uses the provided `Dockerfile` (included in the `offline_bundle`) to create a container image with all dependencies pre-installed.

**Prerequisites for Docker Option:**

*   **Docker Engine Installation:** Docker Engine itself must be installed on the offline Linux machine *before* you can build or run the image. You will need to download the appropriate Docker Engine package(s) (e.g., `.deb` or `.rpm`) and their dependencies for your specific Linux distribution on an *online* machine, transfer them via USB, and install them manually. Refer to the official Docker documentation for offline installation procedures: [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)
*   **Base Image (`ubuntu:20.04`):** The `Dockerfile` relies on the `ubuntu:20.04` base image. This image must also be transferred to the offline machine.
    1.  On an *online* machine with Docker:
        ```bash
        docker pull ubuntu:20.04
        docker save ubuntu:20.04 -o ubuntu2004.tar
        ```
    2.  Transfer `ubuntu2004.tar` via USB to the offline machine.
    3.  On the *offline* machine (after installing Docker Engine):
        ```bash
        docker load -i /path/to/ubuntu2004.tar
        ```

Once Docker Engine is installed and the `ubuntu:20.04` image is loaded on the offline machine, you can proceed:

1.  **Navigate to the bundle directory:**
    ```bash
    cd /path/to/offline_bundle
    ```
    *(The `Dockerfile` is already included in the bundle you transferred).*
2.  **Build the Docker Image:**
    Navigate to the `offline_bundle` directory in the terminal on the Linux machine.
    ```bash
    docker build -t suis_prevalence_analyzer .
    ```
    *(Note: This requires the `ubuntu:20.04` base image to be available locally. If not, it needs to be pulled on an online machine (`docker pull ubuntu:20.04`), saved (`docker save ubuntu:20.04 -o ubuntu2004.tar`), transferred, and loaded (`docker load -i ubuntu2004.tar`) on the offline machine first).*

3.  **Run the Analysis using Docker:**
    ```bash
    docker run --rm suis_prevalence_analyzer
    ```
    This command starts a container from the built image, executes the `run_suis_prevalence.sh` script, and removes the container after completion. Results will be printed to the console. If the script saves files, you might need to mount a volume to access them outside the container (e.g., `docker run --rm -v $(pwd)/results:/app/results suis_prevalence_analyzer`). Adjust the volume mount according to where your script saves output.
