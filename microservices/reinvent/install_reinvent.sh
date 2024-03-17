#!/bin/bash

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup

conda create --name adtools python=2.7 -y
conda install -n adtools conda-forge::xorg-libice -y
conda install -n adtools adfr-suite -c hcc -y

conda create --name reinvent4 python=3.10 -y
#conda install --name reinvent4 numpy swig boost-cpp sphinx sphinx_rtd_theme -y
#conda install --name reinvent4 vina
conda activate reinvent4
pip install -r requirements-linux-64.lock
pip install --no-deps .
