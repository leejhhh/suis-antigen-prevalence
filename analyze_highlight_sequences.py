#!/usr/bin/env python3
"""
S. suis Highlight Sequence Prevalence Analysis
==============================================

이 스크립트는 5개 항원의 highlight 부분 시퀀스를 이용한 prevalence 분석을 수행합니다.
전체 시퀀스 대신 conserved domain이나 중요 영역만을 타겟으로 합니다.

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

def extract_accession(sseqid):
    """Extract genome accession from BLAST subject ID"""
    match = re.search(r'([A-Z]{2}_\d+\.\d+)', sseqid)
    return match.group(1) if match else sseqid

def analyze_highlight_sequences():
    """Highlight 시퀀스들을 분석하는 메인 함수"""
    
    # 설정
    query_fasta = "query_antigens_highlight.fasta"
    output_dir = "suis_highlight_analysis"
    blast_output = f"{output_dir}/blast_results_highlight.tsv"
    db_name = "suis_prevalence_analysis/suis_db"  # 기존 DB 재사용
    total_genomes = 388
    
    # 분석 파라미터
    min_identity = 60.0
    min_coverage = 0.5  # highlight 시퀀스는 더 관대한 coverage 기준 적용
    
    print("🔬 S. suis Highlight Sequence Analysis 시작")
    print("=" * 60)
    
    # 출력 디렉토리 생성
    Path(output_dir).mkdir(exist_ok=True)
    
    # Step 1: Highlight 시퀀스 길이 확인
    print("📏 Highlight 시퀀스 정보 수집...")
    query_lengths = {}
    for seq_record in SeqIO.parse(query_fasta, 'fasta'):
        query_lengths[seq_record.id] = len(seq_record.seq)
        print(f"  {seq_record.id}: {len(seq_record.seq)} aa")
    
    # Step 2: tBLASTn 검색 실행
    print("\n🔍 tBLASTn 검색 실행...")
    cmd = [
        'tblastn',
        '-query', query_fasta,
        '-db', db_name,
        '-evalue', '1e-5',
        '-outfmt', '6',
        '-out', blast_output,
        '-num_threads', '4'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ tBLASTn 실패: {result.stderr}")
        return
    
    print(f"  BLAST 완료: {blast_output}")
    
    # Step 3: 결과 분석
    print("\n📊 BLAST 결과 분석...")
    
    # BLAST 결과 읽기
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore']
    
    try:
        df = pd.read_csv(blast_output, sep='\t', names=cols, header=None)
    except pd.errors.EmptyDataError:
        print("❌ BLAST 결과 파일이 비어있습니다.")
        return
    
    if df.empty:
        print("❌ BLAST hit이 없습니다.")
        return
    
    print(f"  총 BLAST hit 수: {len(df)}")
    
    # 게놈 accession 추가
    df['genome_accession'] = df['sseqid'].apply(extract_accession)
    
    # Coverage 계산
    df['query_length'] = df['qseqid'].map(query_lengths)
    df['coverage'] = df['length'] / df['query_length']
    
    print(f"\n  Raw hit 분포:")
    for antigen in df['qseqid'].unique():
        count = len(df[df['qseqid'] == antigen])
        print(f"    {antigen}: {count} hits")
    
    # 필터링 적용
    filt = df[(df['pident'] >= min_identity) & (df['coverage'] >= min_coverage)]
    print(f"\n  필터링 후 (≥{min_identity}% identity, ≥{min_coverage*100}% coverage): {len(filt)} hits")
    
    # 항원별 분석
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
        
        # 분류
        if prevalence >= 80:
            classification = "High"
        elif prevalence >= 50:
            classification = "Medium"
        else:
            classification = "Low"
        
        # 평가
        if prevalence >= 80:
            assessment = "Excellent highlight region - broad coverage"
        elif prevalence >= 60:
            assessment = "Good highlight region - moderate coverage"
        elif prevalence >= 40:
            assessment = "Moderate highlight region - limited coverage"
        else:
            assessment = "Poor highlight region - very limited"
        
        results.append({
            'antigen': qseqid,
            'highlight_length': query_lengths[qseqid],
            'raw_hits': len(antigen_hits),
            'filtered_hits': len(antigen_filt),
            'hit_genomes': unique_genomes,
            'total_genomes': total_genomes,
            'prevalence_percent': prevalence,
            'max_identity': max_identity,
            'mean_identity': mean_identity,
            'mean_coverage_percent': mean_coverage,
            'classification': classification,
            'assessment': assessment
        })
    
    results_df = pd.DataFrame(results)
    
    # Step 4: 결과 저장
    print("\n💾 결과 저장...")
    
    # TSV 파일 저장
    output_file = f"{output_dir}/highlight_prevalence_stats.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"  결과 저장: {output_file}")
    
    # 요약 보고서 저장
    summary_file = f"{output_dir}/highlight_analysis_summary.txt"
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("S. suis Highlight Sequence Prevalence Analysis\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"분석 파라미터:\n")
        f.write(f"  Identity threshold: ≥{min_identity}%\n")
        f.write(f"  Coverage threshold: ≥{min_coverage*100}%\n")
        f.write(f"  총 genome 수: {total_genomes}\n\n")
        
        f.write("Highlight 시퀀스 결과:\n")
        f.write("-" * 30 + "\n")
        for _, row in results_df.iterrows():
            f.write(f"\n{row['antigen']}:\n")
            f.write(f"  Highlight 길이: {row['highlight_length']} aa\n")
            f.write(f"  Raw hits: {row['raw_hits']}\n")
            f.write(f"  Filtered hits: {row['filtered_hits']}\n")
            f.write(f"  Prevalence: {row['prevalence_percent']:.2f}% ({row['hit_genomes']}/{row['total_genomes']})\n")
            f.write(f"  Max identity: {row['max_identity']:.1f}%\n")
            f.write(f"  Mean identity: {row['mean_identity']:.1f}%\n")
            f.write(f"  Mean coverage: {row['mean_coverage_percent']:.1f}%\n")
            f.write(f"  Classification: {row['classification']}\n")
            f.write(f"  Assessment: {row['assessment']}\n")
    
    print(f"  요약 저장: {summary_file}")
    
    # Step 5: 결과 출력
    print("\n📋 HIGHLIGHT SEQUENCE ANALYSIS 결과:")
    print("-" * 60)
    for _, row in results_df.iterrows():
        antigen_name = row['antigen'].replace('_highlight', '')
        print(f"{antigen_name:12}: {row['prevalence_percent']:6.2f}% ({row['classification']}) - {row['highlight_length']} aa")
    
    print(f"\n✅ Highlight 시퀀스 분석 완료!")
    print(f"📊 결과 디렉토리: {output_dir}/")
    
    return results_df

if __name__ == "__main__":
    analyze_highlight_sequences() 