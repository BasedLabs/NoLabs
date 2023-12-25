import glob
import os
import shutil
import subprocess
import traceback
from typing import List, Dict

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
rfdiffusion_dir = os.path.join(dname, 'RFdiffusion')
os.chdir(rfdiffusion_dir)
input_pdbs_dir = os.path.join(dname, 'input_pdbs')
output_files_dir = os.path.join(dname, 'output_pdbs')

pdb_tmp_file_path = os.path.join(input_pdbs_dir, 'tmp.pdb')


def ensure_dirs_exist():
    if not os.path.exists(input_pdbs_dir):
        os.mkdir(input_pdbs_dir)

    if not os.path.exists(output_files_dir):
        os.mkdir(output_files_dir)


def cleanup_folders():
    shutil.rmtree(input_pdbs_dir)
    shutil.rmtree(output_files_dir)


def pipeline(pdb_content: str = None,
             contig: str = '50',
             timesteps: int = None,
             hotspots: str = None,
             number_of_designs: int = 1) -> Dict:
    """
                    Run raw inference on rfdiffusion
                    :param pdb_content: pdb content (optional). If not specified design an unconditional monomer
                    :param contig: contigs (around what and how to design a protein). See https://github.com/RosettaCommons/RFdiffusion
                    Default None
                    :param timesteps: default 50 ( Desired iterations to generate structure. )
                    :param hotspots: A30, A33, A34
                    :param number_of_designs: Number of designs to generate
                    The model optionally readily learns that it should be making an interface which involving these hotspot residues. Input is ChainResidueNumber: A100 for residue 100 on chain A.
                    :return:
                    """
    try:
        ensure_dirs_exist()
        inference_path = os.path.join(rfdiffusion_dir, 'scripts', 'run_inference.py')
        program = ['python3.9', inference_path, 'inference.model_directory_path=/app/RFdiffusion/models',
                   f'inference.num_designs={number_of_designs}']

        program.append(f'inference.output_prefix={output_files_dir}/result')

        if pdb_content:
            with open(pdb_tmp_file_path, 'w') as f:
                f.write(pdb_content)
            program.append(f'inference.input_pdb={pdb_tmp_file_path}')

        program.append(f'contigmap.contigs=[{contig}]')

        program.append('denoiser.noise_scale_ca=0')
        program.append('denoiser.noise_scale_frame=0')

        if timesteps:
            program.append(f'diffuser.T={timesteps}')

        if hotspots:
            program.append(f'ppi.hotspot_res=[{hotspots}]')

        res = subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        print(res.stdout.decode('utf-8'))
        print(res.stderr.decode('utf-8'))

        pdbs = []
        files = glob.glob(os.path.join(output_files_dir, '*.pdb'))
        for file_name in files:
            with open(file_name, 'r') as f:
                pdbs.append(f.read())
        return {'pdbs': pdbs, 'errors': []}
    except ValueError:
        return {'pdbs': [], 'errors': ['Contig is incorrect']}
    except Exception:
        return {'pdbs': [], 'errors': ['Unknown error']}
    # TODO add exception handling for incorrect hotspots
    finally:
        cleanup_folders()
