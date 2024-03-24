#!/bin/bash

# make the script stop when error (non-true exit code) occurs
set -e

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

SCRIPT=`realpath -s $0`
export PIPEDIR=`dirname $SCRIPT`
HHDB="$PIPEDIR/pdb100_2021Mar03/pdb100_2021Mar03"

conda activate RF2

echo '$1 - msa file, $2 prefix to output pdb'

python $PIPEDIR/network/predict.py -inputs $1 -prefix $2 -model $PIPEDIR/network/weights/RF2_apr23.pt -db $PIPEDIR/pdb100_2021Mar03/pdb100_2021Mar03 -symm C1