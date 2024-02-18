import os
import json
import tempfile

# update these paths to reflect your system's configuration
dockstream_path = os.path.expanduser("~/Desktop/ProjectData/DockStream")
dockstream_env = os.path.expanduser("~/miniconda3/envs/DockStream")
vina_binary_location = os.path.expanduser(
    "~/Desktop/ProjectData/foreign/AutoDockVina/autodock_vina_1_1_2_linux_x86/bin")

# no changes are necessary beyond this point
# ---------
# get the notebook's root path
try:
    ipynb_path
except NameError:
    ipynb_path = os.getcwd()

# generate the paths to the entry points
target_preparator = dockstream_path + "/target_preparator.py"
docker = dockstream_path + "/docker.py"

# generate a folder to store the results
output_dir = os.path.expanduser("~/Desktop/AutoDock_Vina_demo")
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

# generate the paths to the files shipped with this implementation
apo_1UYD_path = ipynb_path + "/../data/1UYD/1UYD_apo.pdb"
reference_ligand_path = ipynb_path + "/../data/1UYD/PU8.pdb"
smiles_path = ipynb_path + "/../data/1UYD/ligands_smiles.txt"

# generate output paths for the configuration file, the "fixed" PDB file and the "Gold" receptor
target_prep_path = output_dir + "/ADV_target_prep.json"
fixed_pdb_path = output_dir + "/ADV_fixed_target.pdb"
adv_receptor_path = output_dir + "/ADV_receptor.pdbqt"
log_file_target_prep = output_dir + "/ADV_target_prep.log"
log_file_docking = output_dir + "/ADV_docking.log"

# generate output paths for the configuration file, embedded ligands, the docked ligands and the scores
docking_path = output_dir + "/ADV_docking.json"
ligands_conformers_path = output_dir + "/ADV_embedded_ligands.sdf"
ligands_docked_path = output_dir + "/ADV_ligands_docked.sdf"
ligands_scores_path = output_dir + "/ADV_scores.csv"

# -------------------------------------------------------------------------------------------------------

# specify the target preparation JSON file as a dictionary and write it out
tp_dict = {
    "target_preparation":
        {
            "header": {  # general settings
                "logging": {  # logging settings (e.g. which file to write to)
                    "logfile": log_file_target_prep
                }
            },
            "input_path": apo_1UYD_path,  # this should be an absolute path
            "fixer": {  # based on "PDBFixer"; tries to fix common problems with PDB files
                "enabled": True,
                "standardize": True,  # enables standardization of residues
                "remove_heterogens": True,  # remove hetero-entries
                "fix_missing_heavy_atoms": True,  # if possible, fix missing heavy atoms
                "fix_missing_hydrogens": True,  # add hydrogens, which are usually not present in PDB files
                "fix_missing_loops": False,  # add missing loops; CAUTION: the result is usually not sufficient
                "add_water_box": False,  # if you want to put the receptor into a box of water molecules
                "fixed_pdb_path": fixed_pdb_path  # if specified and not "None", the fixed PDB file will be stored here
            },
            "runs": [  # "runs" holds a list of backend runs; at least one is required
                {
                    "backend": "AutoDockVina",  # one of the backends supported ("AutoDockVina", "OpenEye", ...)
                    "output": {
                        "receptor_path": adv_receptor_path  # the generated receptor file will be saved to this location
                    },
                    "parameters": {
                        "pH": 7.4,  # sets the protonation states (NOT used in Vina)
                        "extract_box": {  # in order to extract the coordinates of the pocket (see text)
                            "reference_ligand_path": reference_ligand_path,  # path to the reference ligand
                            "reference_ligand_format": "PDB"  # format of the reference ligand
                        }
                    }}]}}

with open(target_prep_path, 'w') as f:
    json.dump(tp_dict, f, indent="    ")

# -------- DOCKING --------------
# load the smiles (just for illustrative purposes)
# here, 15 moleucles will be used
with open(smiles_path, 'r') as f:
    smiles = [smile.strip() for smile in f.readlines()]
print(smiles)

# specify the embedding and docking JSON file as a dictionary and write it out
ed_dict = {
    "docking": {
        "header": {  # general settings
            "logging": {  # logging settings (e.g. which file to write to)
                "logfile": log_file_docking
            }
        },
        "ligand_preparation": {  # the ligand preparation part, defines how to build the pool
            "embedding_pools": [
                {
                    "pool_id": "Corina_pool",  # here, we only have one pool
                    "type": "Corina",
                    "parameters": {
                        "prefix_execution": "module load corina"
                        # only required, if a module needs to be loaded to execute "Corina"
                    },
                    "input": {
                        "standardize_smiles": False,
                        "type": "smi",
                        "input_path": smiles_path
                    },
                    "output": {  # the conformers can be written to a file, but "output" is
                        # not required as the ligands are forwarded internally
                        "conformer_path": ligands_conformers_path,
                        "format": "sdf"
                    }
                }
            ]
        },
        "docking_runs": [
            {
                "backend": "AutoDockVina",
                "run_id": "AutoDockVina",
                "input_pools": ["Corina_pool"],
                "parameters": {
                    "binary_location": vina_binary_location,  # absolute path to the folder, where the "vina" binary
                    # can be found
                    "parallelization": {
                        "number_cores": 4
                    },
                    "seed": 42,  # use this "seed" to generate reproducible results; if
                    # varied, slightly different results will be produced
                    "receptor_pdbqt_path": [adv_receptor_path],  # paths to the receptor files
                    "number_poses": 2,  # number of poses to be generated
                    "search_space": {  # search space (cavity definition); see text
                        "--center_x": 3.3,
                        "--center_y": 11.5,
                        "--center_z": 24.8,
                        "--size_x": 15,
                        "--size_y": 10,
                        "--size_z": 10
                    }
                },
                "output": {
                    "poses": {"poses_path": ligands_docked_path},
                    "scores": {"scores_path": ligands_scores_path}
                }}]}}

with open(docking_path, 'w') as f:
    json.dump(ed_dict, f, indent=2)

# print out path to generated JSON
print(docking_path)

# ------ Docking backend ---------

configuration = {'scoring_function': {}}  # type: ignore
configuration['scoring_function']['parameters'] = {
    "component_type": "dockstream",
    "name": "dockstream",
    "weight": 1,
    "specific_parameters": {
        "transformation": {
            "transformation_type": "reverse_sigmoid",
            "low": -12,
            "high": -8,
            "k": 0.25
        },
        "configuration_path": "<absolute_path_to_DockStream_configuration>/docking.json",
        "docker_script_path": "<absolute_path_to_DockStream_source>/docker.py",
        "environment_path": "<absolute_path_to_miniconda_installation>/envs/DockStream/bin/python"
    }
}
