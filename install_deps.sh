#!/bin/bash
pip install --upgrade pip
pip install torch==1.13.1 torchvision torchaudio
pip install torch-geometric
pip install --no-index pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.0+cpu.html
pip install -r requirements.txt