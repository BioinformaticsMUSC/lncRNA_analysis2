import os
import re
import glob
import argparse
import pandas as pd
from collections import defaultdict

# regex pattern to detect _R1/_1/_R2/_2 just before .fastq.gz
PAIR_PATTERN = re.compile(r'(.*?)([_\.](R?1|R?2))\.f(ast)?q\.gz$', re.IGNORECASE)

def find_fastq_files(base_dir):
    fastq_files = glob.glob(os.path.join(base_dir, '**', '*.fastq.gz'), recursive=True)
    sample_map = defaultdict(dict)

    for fq in fastq_files:
        fname = os.path.basename(fq)
        match = PAIR_PATTERN.match(fname)
        
        if match:
            sample_id, _, read_number = match.groups()
            read_number = read_number.lower()
            if '1' in read_number:
                sample_map[sample_id]['fq1'] = fq
            elif '2' in read_number:
                sample_map[sample_id]['fq2'] = fq
        else:
            # No identifiable _1/_2, assume single-end
            sample_id = os.path.splitext(os.path.splitext(fname)[0])[0]  # remove .fastq.gz
            sample_map[sample_id]['fq1'] = fq

    return sample_map

def write_samples_tsv(sample_map, output_file):
    rows = []
    for sample, reads in sorted(sample_map.items()):
        rows.append({
            'sample': sample,
            'fq1': reads.get('fq1', ''),
            'fq2': reads.get('fq2', '')
        })
    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Robustly generate samples.tsv from FASTQ files')
    parser.add_argument('indir', help='Input directory containing FASTQ files (recursive)')
    parser.add_argument('-o', '--outfile', default='samples.tsv', help='Output TSV file (default: samples.tsv)')
    args = parser.parse_args()

    sample_map = find_fastq_files(args.indir)
    write_samples_tsv(sample_map, args.outfile)
    print(f"✅ Done. Found {len(sample_map)} samples. Output written to {args.outfile}")
