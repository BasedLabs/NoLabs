import csv
import sys

samping_direct = sys.argv[1] # sampling_direct.csv
scoring_input = sys.argv[2] # scoring_input.smi

with open(samping_direct, 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

smiles = [r['SMILES'] for r in rows]
with open(scoring_input, 'w') as f:
    for smi in smiles:
        f.write(smi + '\n')