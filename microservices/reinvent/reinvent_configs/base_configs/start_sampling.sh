#!/bin/bash

echo "Usage: ./start_sampling.sh Sampling_file_path Scoring_file_path log_file_output error_file_output"

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
conda activate reinvent4

BASEDIR=$(dirname "$0")

reinvent $1 1>> $3 2>> $4

python $BASEDIR/extract_smi.py "$BASEDIR/sampling_direct.csv" "$BASEDIR/scoring_input.smi"

reinvent $2 1>> $3 2>> $4