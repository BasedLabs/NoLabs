#!/bin/bash

echo "Usage: ./start_reinforcement_learning.sh RL_file_path log_file_output error_file_output"

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
conda activate reinvent4

BASEDIR=$(dirname "$0")

reinvent $1 1> $2 2> $3

