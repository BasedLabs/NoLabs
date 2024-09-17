#!/bin/bash

echo "Usage: ./prepare_receptor.sh input-pdb-file output-pdbqt-file"

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
conda activate adtools
prepare_receptor -r $1 -o $2
