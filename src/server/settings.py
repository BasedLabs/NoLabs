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
CONFORMATIONS_EXPERIMENTS_DIR = EXPERIMENTS_DIR + "/conformations"
FASTA_API = 'http://207.246.89.242:8000/generate-msa/'
