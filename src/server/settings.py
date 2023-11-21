import os

is_test = True
use_gpu = False
host = '127.0.0.1'
port = 5000

dirname = os.path.dirname
# Base directory for storing experiment results
EXPERIMENTS_DIR = dirname(dirname(os.path.abspath(__file__)) + "/experiments/")
PROTEIN_EXPERIMENTS_DIR = EXPERIMENTS_DIR + "/proteins"
DTI_EXPERIMENTS_DIR = EXPERIMENTS_DIR + "/drug_discovery"
FASTA_API = 'http://143.198.141.110:8000/generate-msa/'
CONFORMATIONS_EXPERIMENTS_DIR = EXPERIMENTS_DIR + "/conformations"