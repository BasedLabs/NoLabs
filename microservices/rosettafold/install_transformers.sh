#!/bin/bash

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup

conda activate RF2
cd SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
pip install torchdata