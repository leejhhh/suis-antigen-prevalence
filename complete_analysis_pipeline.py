#!/usr/bin/env python3
"""
S. suis 5-Antigen Prevalence Analysis - Complete Pipeline
=======================================================

This script reproduces the complete analysis pipeline used for the 
S. suis antigen prevalence study reported in the accompanying manuscript.

Requirements:
- BLAST+ (version 2.15.0 or later)
- Python 3.x with pandas, biopython, numpy
- Input files: query_antigens.fasta, suis_selected/ directory

Usage:
    python complete_analysis_pipeline.py

Author: [Principal Investigator]
Date: May 26, 2025
"""

import os
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import re
from pathlib import Path

class SsuisAntiGenAnalyzer:
    def __init__(self, config=None):
        """Initialize analyzer with configuration parameters"""
        self.config = config or {
            'query_fasta': 'query_antigens.fasta',
            'genome_dir': 'suis_selected',
            'output_dir': 'suis_prevalence_analysis', 
            'evalue': '1e-5',
            'min_identity': 60.0,
            'min_coverage': 0.8,
            'threads': 4
        }
        
        # Create output directory
        Path(self.config['output_dir']).mkdir(exist_ok=True)
        
        # Expected antigen information (from analysis)
        self.antigen_info = {
            'HP0197|WP_277937340.1': {'length': 671, 'function': 'Hypothetical protein'},
            'Fnb|WP_014636551.1': {'length': 577, 'function': 'Fibronectin-binding protein'},
            'SAO|WP_211840080.1': {'length': 407, 'function': 'Surface antigen one'},
            'C5a|WP_240208248.1': {'length': 492, 'function': 'C5a peptidase'},
            'Suilysin|AIG43067.1': {'length': 475, 'function': 'Cytolysin'}
        }
    
    def extract_accession(self, sseqid):
        """Extract genome accession from BLAST subject ID"""
        match = re.search(r'([A-Z]{2}_\d+\.\d+)', sseqid)
        return match.group(1) if match else sseqid
    
    def validate_inputs(self):
        """Validate required input files and tools"""
        print("🔍 Validating inputs...")
        
        # Check query FASTA
        if not os.path.exists(self.config['query_fasta']):
            raise FileNotFoundError(f"Query FASTA not found: {self.config['query_fasta']}")
        
        # Check genome directory
        genome_dir = Path(self.config['genome_dir'])
        if not genome_dir.exists():
            raise FileNotFoundError(f"Genome directory not found: {genome_dir}")
        
        # Count genome files
        genome_files = list(genome_dir.glob('*.fna'))
        print(f"  Found {len(genome_files)} genome files")
        
        # Check BLAST tools
        try:
            subprocess.run(['tblastn', '-version'], capture_output=True, check=True)
            subprocess.run(['makeblastdb', '-version'], capture_output=True, check=True)
            print("  BLAST+ tools: OK")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError("BLAST+ tools not found in PATH")
        
        return len(genome_files)
    
    def merge_genomes(self):
        """Merge all genome files into single FASTA"""
        print("🧬 Merging genome files...")
        
        merged_file = Path(self.config['output_dir']) / 'all_suis_genomes.fna'
        genome_dir = Path(self.config['genome_dir'])
        
        # Use system command for efficiency (works on both Windows/Linux)
        if os.name == 'nt':  # Windows
            cmd = f'type "{genome_dir}\\*.fna" > "{merged_file}"'
            subprocess.run(cmd, shell=True, check=True)
        else:  # Unix/Linux
            cmd = f'cat "{genome_dir}"/*.fna > "{merged_file}"'
            subprocess.run(cmd, shell=True, check=True)
        
        print(f"  Merged file created: {merged_file}")
        return merged_file
    
    def create_blast_database(self, merged_fasta):
        """Create BLAST nucleotide database"""
        print("🗃️ Creating BLAST database...")
        
        db_name = Path(self.config['output_dir']) / 'suis_db'
        
        cmd = [
            'makeblastdb',
            '-in', str(merged_fasta),
            '-dbtype', 'nucl',
            '-out', str(db_name),
            '-parse_seqids',
            '-title', 'S.suis_DB'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"makeblastdb failed: {result.stderr}")
        
        print(f"  Database created: {db_name}")
        return db_name
    
    def run_tblastn_search(self, db_name):
        """Run tBLASTn search against database"""
        print("🔬 Running tBLASTn search...")
        
        blast_output = Path(self.config['output_dir']) / 'blast_results.tsv'
        
        cmd = [
            'tblastn',
            '-query', self.config['query_fasta'],
            '-db', str(db_name),
            '-evalue', self.config['evalue'],
            '-outfmt', '6',  # Tabular format
            '-out', str(blast_output),
            '-num_threads', str(self.config['threads'])
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"tblastn failed: {result.stderr}")
        
        print(f"  BLAST completed: {blast_output}")
        return blast_output
    
    def load_query_lengths(self):
        """Load query protein lengths from FASTA file"""
        print("📏 Loading query protein lengths...")
        
        query_lengths = {}
        for seq_record in SeqIO.parse(self.config['query_fasta'], 'fasta'):
            query_lengths[seq_record.id] = len(seq_record.seq)
            print(f"  {seq_record.id}: {len(seq_record.seq)} aa")
        
        return query_lengths
    
    def analyze_blast_results(self, blast_file, total_genomes, query_lengths):
        """Analyze BLAST results and calculate prevalence"""
        print("📊 Analyzing BLAST results...")
        
        # Load BLAST results
        cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
                'qstart','qend','sstart','send','evalue','bitscore']
        
        try:
            df = pd.read_csv(blast_file, sep='\t', names=cols, header=None)
        except pd.errors.EmptyDataError:
            print("  Warning: BLAST results file is empty")
            return pd.DataFrame()
        
        if df.empty:
            print("  Warning: No BLAST hits found")
            return pd.DataFrame()
        
        print(f"  Total BLAST hits: {len(df)}")
        
        # Add genome accession
        df['genome_accession'] = df['sseqid'].apply(self.extract_accession)
        
        # Calculate coverage
        df['query_length'] = df['qseqid'].map(query_lengths)
        df['coverage'] = df['length'] / df['query_length']
        
        print(f"  Raw hit distribution:")
        for antigen in df['qseqid'].unique():
            count = len(df[df['qseqid'] == antigen])
            print(f"    {antigen}: {count} hits")
        
        # Apply filters
        min_identity = self.config['min_identity']
        min_coverage = self.config['min_coverage']
        
        filt = df[(df['pident'] >= min_identity) & (df['coverage'] >= min_coverage)]
        print(f"  Hits after filtering (≥{min_identity}% identity, ≥{min_coverage*100}% coverage): {len(filt)}")
        
        # Analyze by antigen
        results = []
        
        for qseqid in query_lengths.keys():
            antigen_hits = df[df['qseqid'] == qseqid]
            antigen_filt = filt[filt['qseqid'] == qseqid]
            
            if len(antigen_filt) > 0:
                unique_genomes = antigen_filt['genome_accession'].nunique()
                prevalence = (unique_genomes / total_genomes) * 100
                max_identity = antigen_filt['pident'].max()
                mean_identity = antigen_filt['pident'].mean()
                mean_coverage = antigen_filt['coverage'].mean() * 100
            else:
                unique_genomes = 0
                prevalence = 0.0
                max_identity = 0.0
                mean_identity = 0.0
                mean_coverage = 0.0
            
            results.append({
                'antigen': qseqid,
                'protein_length': query_lengths[qseqid],
                'raw_hits': len(antigen_hits),
                'filtered_hits': len(antigen_filt),
                'hit_genomes': unique_genomes,
                'total_genomes': total_genomes,
                'prevalence_percent': prevalence,
                'max_identity': max_identity,
                'mean_identity': mean_identity,
                'mean_coverage_percent': mean_coverage,
                'classification': self._classify_prevalence(prevalence),
                'assessment': self._assess_vaccine_potential(qseqid, prevalence, len(antigen_hits))
            })
        
        return pd.DataFrame(results)
    
    def _classify_prevalence(self, prevalence):
        """Classify prevalence level"""
        if prevalence >= 80:
            return "High"
        elif prevalence >= 50:
            return "Medium"
        else:
            return "Low"
    
    def _assess_vaccine_potential(self, antigen, prevalence, raw_hits):
        """Assess vaccine potential based on prevalence and other factors"""
        if 'HP0197' in antigen and prevalence == 0 and raw_hits > 0:
            return "Coverage threshold issue - requires re-analysis"
        elif prevalence >= 80:
            return "Excellent vaccine candidate - broad coverage"
        elif prevalence >= 60:
            return "Moderate vaccine candidate - regional variations possible"
        elif prevalence >= 40:
            return "Moderate coverage - highly conserved when present"
        else:
            return "Poor vaccine candidate - limited distribution"
    
    def save_results(self, results_df):
        """Save analysis results"""
        print("💾 Saving results...")
        
        output_file = Path(self.config['output_dir']) / 'detailed_antigen_stats.tsv'
        results_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Results saved: {output_file}")
        
        # Also save formatted summary
        summary_file = Path(self.config['output_dir']) / 'analysis_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("S. suis 5-Antigen Prevalence Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Analysis parameters:\n")
            f.write(f"  Identity threshold: ≥{self.config['min_identity']}%\n")
            f.write(f"  Coverage threshold: ≥{self.config['min_coverage']*100}%\n")
            f.write(f"  E-value threshold: {self.config['evalue']}\n\n")
            
            for _, row in results_df.iterrows():
                f.write(f"{row['antigen']}:\n")
                f.write(f"  Prevalence: {row['prevalence_percent']:.2f}% ({row['hit_genomes']}/{row['total_genomes']})\n")
                f.write(f"  Classification: {row['classification']}\n")
                f.write(f"  Assessment: {row['assessment']}\n\n")
        
        print(f"  Summary saved: {summary_file}")
        return output_file
    
    def run_complete_analysis(self):
        """Run the complete analysis pipeline"""
        print("🚀 Starting S. suis antigen prevalence analysis...")
        print("=" * 60)
        
        try:
            # Step 1: Validate inputs
            total_genomes = self.validate_inputs()
            
            # Step 2: Merge genomes
            merged_fasta = self.merge_genomes()
            
            # Step 3: Create BLAST database
            db_name = self.create_blast_database(merged_fasta)
            
            # Step 4: Run BLAST search
            blast_output = self.run_tblastn_search(db_name)
            
            # Step 5: Load query lengths
            query_lengths = self.load_query_lengths()
            
            # Step 6: Analyze results
            results_df = self.analyze_blast_results(blast_output, total_genomes, query_lengths)
            
            # Step 7: Save results
            output_file = self.save_results(results_df)
            
            print("✅ Analysis completed successfully!")
            print(f"📊 Results available in: {self.config['output_dir']}")
            
            return results_df
            
        except Exception as e:
            print(f"❌ Analysis failed: {e}")
            raise

def main():
    """Main function to run the analysis"""
    
    # Configuration (modify as needed)
    config = {
        'query_fasta': 'query_antigens.fasta',
        'genome_dir': 'suis_selected',
        'output_dir': 'suis_prevalence_analysis',
        'evalue': '1e-5',
        'min_identity': 60.0,
        'min_coverage': 0.8,
        'threads': 4
    }
    
    # Run analysis
    analyzer = SsuisAntiGenAnalyzer(config)
    results = analyzer.run_complete_analysis()
    
    # Display summary
    print("\n📋 FINAL SUMMARY:")
    print("-" * 40)
    for _, row in results.iterrows():
        antigen_name = row['antigen'].split('|')[0]
        print(f"{antigen_name:12}: {row['prevalence_percent']:6.2f}% ({row['classification']})")

if __name__ == "__main__":
    main() 