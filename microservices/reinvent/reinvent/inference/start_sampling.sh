#!/bin/bash

echo "Usage: ./start_sampling.sh sampling_toml_path"

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
conda activate reinvent4
reinvent $1
