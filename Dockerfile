# Use a base image with common tools
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install basic dependencies
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    unzip \
    procps \
    jq \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- Install Miniconda ---
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -a -y && \
    ln -s $CONDA_DIR/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    find $CONDA_DIR -follow -type f -name '*.js.map' -delete

# --- Install BLAST+ ---
# Download specific version for reproducibility and offline use preparation
ARG BLAST_VERSION=2.15.0
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    tar zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    mv ncbi-blast-${BLAST_VERSION}+/bin/* /usr/local/bin/ && \
    rm -rf ncbi-blast-${BLAST_VERSION}+*

# --- Install NCBI Datasets CLI ---
# Download latest stable version from NCBI FTP
# Install curl temporarily, download, make executable, then remove curl
RUN apt-get update && apt-get install -y --no-install-recommends curl && \
    curl -o /usr/local/bin/datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' && \
    chmod +x /usr/local/bin/datasets && \
    apt-get purge -y --auto-remove curl && \
    rm -rf /var/lib/apt/lists/*

# --- Create Conda Environment ---
# Create environment with specified Python packages
RUN conda create -y -n suis_env python=3.10 pandas biopython && \
    conda clean -afy

# Activate conda environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate suis_env" >> ~/.bashrc
ENV PATH="/opt/conda/envs/suis_env/bin:$PATH"

# --- Setup Application ---
WORKDIR /app

# Copy pipeline scripts
COPY run_suis_prevalence.sh .
COPY parse_prevalence.py .

# Ensure scripts are executable
RUN chmod +x run_suis_prevalence.sh

# Create output directory structure expected by the script
RUN mkdir -p suis_prevalence_analysis/genomes

# Set entrypoint or default command (optional)
# The user will typically run the script manually after providing the query fasta
# ENTRYPOINT ["./run_suis_prevalence.sh"]

# Add a note about providing the query file
CMD echo "Container ready. \n1. Copy your query protein FASTA file to the container's /app directory as 'query_protein.fasta'. \n   Example: docker cp query_protein.fasta <container_id>:/app/query_protein.fasta \n2. Execute the pipeline: docker exec -it <container_id> ./run_suis_prevalence.sh"
