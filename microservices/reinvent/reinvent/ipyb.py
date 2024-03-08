# load dependencies
import os
import re
import json
import tempfile

# --------- change these path variables as required
reinvent_dir = os.path.expanduser("/app/REINVENT4")
reinvent_env = os.path.expanduser("/var/conda/envs/reinvent4")

# DockStream variables
dockstream_dir = os.path.expanduser("/app/DockStream")
dockstream_env = os.path.expanduser("/var/conda/envs/DockStream")
# generate the path to the DockStream entry points
docker_path = os.path.join(dockstream_dir, "docker.py")

output_dir = os.path.expanduser("/app/demo")

# --------- do not change
# get the notebook's root path
try: ipynb_path
except NameError: ipynb_path = os.getcwd()

# if required, generate the folder to store the results
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

# Glide docking variables
grid_file_path = os.path.expanduser("/app/notebooks/data/DockStream/1UYD_grid.zip")
output_ligands_docked_poses_path = os.path.expanduser("/app/demo/docked_poses")
output_ligands_docking_scores_path = os.path.expanduser("/app/demo/docking_scores")

try:
    os.mkdir(output_ligands_docked_poses_path)
except FileExistsError:
    pass

try:
    os.mkdir(output_ligands_docking_scores_path)
except FileExistsError:
    pass

docking_configuration_path = os.path.join(output_dir, "Glide_DockStream_Conf.json")