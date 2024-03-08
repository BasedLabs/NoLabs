#!/bin/bash

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup

conda create --name reinvent4 python=3.10 -y
conda activate reinvent4
pip install -r requirements-linux-64.lock
pip install --no-deps .
pip install adfr-suite