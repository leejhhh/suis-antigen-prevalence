import subprocess
from pathlib import Path

def test_toy_blast():
    root = Path(__file__).resolve().parents[1]
    blast_file = root / 'sample_data' / 'toy_blast.tsv'
    query = root / 'query_antigens.fasta'
    out_tsv = root / 'sample_data' / 'toy_stats.tsv'

    cmd = [
        'python', str(root / 'parse_prevalence.py'),
        '-i', str(blast_file),
        '-t', '3',
        '-q', str(query),
        '-o', str(out_tsv)
    ]
    subprocess.check_call(cmd)

    assert out_tsv.exists(), "Output TSV not created"
    with out_tsv.open() as f:
        lines = f.readlines()
    # header + 3 lines expected
    assert len(lines) == 4, f"Expected 4 lines, got {len(lines)}"
