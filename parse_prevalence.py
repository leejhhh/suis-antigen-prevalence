#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from Bio import SeqIO

# Helper function to safely extract accession
def extract_accession(sseqid):
    match = re.search(r'([A-Z]{2}_\d+\.\d+)', sseqid)
    return match.group(1) if match else sseqid # Return original sseqid if no match

def parse_blast_output(blast_file, total_genomes, query_fasta,
                       min_identity, min_coverage):
    """
    Parses BLAST tabular output (format 6), filters by identity/coverage,
    and calculates prevalence statistics.

    Args:
        blast_file (str): Path to the BLAST output file (fmt 6).
        total_genomes (int): Total number of genomes from metadata.
        query_fasta (str): Path to the query protein FASTA file.
        min_identity (float): Minimum percent identity threshold.
        min_coverage (float): Minimum query coverage threshold (fraction, e.g., 0.8).

    Returns:
        tuple: (prevalence_percentage, hit_stats_df)
               - prevalence_percentage (float): Percentage of genomes with hits passing filters.
               - hit_stats_df (pd.DataFrame): DataFrame with stats per hit genome.
    """
    try:
        # Read BLAST results
        cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
                'qstart','qend','sstart','send','evalue','bitscore']
        df = pd.read_csv(blast_file, sep='\t', names=cols, header=None)

        if df.empty:
            print("No hits found in BLAST results.")
            return 0.0, pd.DataFrame()

        # Calculate Query length from the first sequence in the FASTA
        try:
            query_seq = next(SeqIO.parse(query_fasta, 'fasta'))
            qlen = len(query_seq.seq)
            print(f"Query sequence ID: {query_seq.id}, Length: {qlen}")
        except FileNotFoundError:
            print(f"Error: Query FASTA file not found at {query_fasta}")
            exit(1)
        except StopIteration:
             print(f"Error: No sequences found in query FASTA file {query_fasta}")
             exit(1)
        except Exception as e:
            print(f"Error reading query FASTA file {query_fasta}: {e}")
            exit(1)

        if qlen == 0:
            print("Error: Query sequence length is 0.")
            exit(1)

        # Add genome accession column
        df['genome_accession'] = df['sseqid'].apply(extract_accession)
        # Calculate coverage (alignment length / query length)
        # Note: 'length' in BLAST fmt 6 is the alignment length
        df['coverage'] = df['length'] / qlen

        # Apply identity and coverage filters
        print(f"Applying filters: Identity >= {min_identity}%, Coverage >= {min_coverage*100:.1f}%")
        filt = df[(df['pident'] >= min_identity) & (df['coverage'] >= min_coverage)].copy() # Use .copy() to avoid SettingWithCopyWarning

        if filt.empty:
            print("No hits passed the identity/coverage filters.")
            return 0.0, pd.DataFrame()

        # Aggregate stats per genome: count hits, max identity, mean coverage
        print("Aggregating statistics per genome...")
        grp = filt.groupby('genome_accession')
        stats = pd.DataFrame({
            'hit_count': grp.size(),
            'max_identity': grp['pident'].max(),
            'mean_coverage': grp['coverage'].mean() # Mean coverage of hits passing filters
        }).reset_index()

        # Calculate prevalence percentage
        num_hit_genomes = stats.shape[0]
        prevalence_percentage = (num_hit_genomes / total_genomes) * 100 if total_genomes > 0 else 0

        # Format output DataFrame
        stats['max_identity'] = stats['max_identity'].round(1)
        stats['mean_coverage'] = (stats['mean_coverage'] * 100).round(1) # Convert coverage back to percentage for output

        print(f"\n--- Results Summary ---")
        print(f"Total genomes analyzed (from metadata): {total_genomes}")
        print(f"Genomes with hits passing filters: {num_hit_genomes}")
        print(f"Prevalence (filtered): {prevalence_percentage:.2f}%")
        print(f"-----------------------\n")

        return prevalence_percentage, stats

    except FileNotFoundError:
        print(f"Error: BLAST output file not found at {blast_file}")
        exit(1)
    except Exception as e:
        print(f"An error occurred during BLAST parsing: {e}")
        exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse BLAST output (fmt 6), filter by identity/coverage, and calculate protein prevalence statistics."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the BLAST output file (tabular format 6).")
    parser.add_argument("-t", "--total_genomes", required=True, type=int, help="Total number of genomes (from metadata JSON).")
    parser.add_argument("-q", "--query_fasta", required=True, help="Path to the query protein FASTA file (used for length calculation).")
    parser.add_argument("--min_identity", type=float, default=70.0, help="Minimum percent identity threshold (default: 70.0).")
    parser.add_argument("--min_coverage", type=float, default=0.8, help="Minimum query coverage threshold (fraction, e.g., 0.8 for 80%, default: 0.8).")
    parser.add_argument("-o", "--output", default="genomes_with_hit_stats.tsv", help="Path to save the filtered hit statistics (TSV format). Default: genomes_with_hit_stats.tsv")

    args = parser.parse_args()

    # Validate coverage input
    if not 0.0 <= args.min_coverage <= 1.0:
        print("Error: --min_coverage must be between 0.0 and 1.0.")
        exit(1)

    prevalence, df_stats = parse_blast_output(
        args.input, args.total_genomes, args.query_fasta,
        args.min_identity, args.min_coverage
    )

    # Save the statistics DataFrame
    if not df_stats.empty:
        # Define output columns explicitly for order
        output_cols = ['genome_accession', 'hit_count', 'max_identity', 'mean_coverage']
        df_stats[output_cols].to_csv(args.output, sep='\t', index=False, float_format='%.1f')
        print(f"Filtered hit statistics saved to: {args.output}")
    else:
        # Create an empty file with header if no hits passed filters
        with open(args.output, 'w') as f:
            f.write("genome_accession\thit_count\tmax_identity\tmean_coverage\n")
        print(f"No hits passed filters. Empty file with header created: {args.output}")
