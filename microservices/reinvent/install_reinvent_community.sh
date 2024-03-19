#!/bin/bash

set -e

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup

conda activate ReinventCommunity
wget https://bootstrap.pypa.io/get-pip.py -O get-pip.py
python get-pip.py --force-reinstall
pip install pyqt5-sip==4.19.18
pip install pyqtwebengine==5.12.1