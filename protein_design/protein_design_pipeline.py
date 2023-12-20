import glob
import os
import subprocess
import traceback

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
rfdiffusion_dir = os.path.join(dname, 'RFdiffusion')
os.chdir(rfdiffusion_dir)
input_pdbs_dir = os.path.join(dname, 'input_pdbs')
output_files_dir = os.path.join(dname, 'output_pdbs')

pdb_tmp_file = 'tmp.pdb'
output_tmp_file = 'output.pdb'


def cleanup_folders():
    files = glob.glob(os.path.join(input_pdbs_dir, '*'))
    for f in files:
        os.remove(f)
    files = glob.glob(os.path.join(output_files_dir, '*'))
    for f in files:
        os.remove(f)


def pipeline(pdb_content: str = None,
             contig: str = '50',
             symmetry: str = None,
             timesteps: int = None,
             hotspots: str = None):
    """
                    Run raw inference on rfdiffusion
                    :param pdb_content: pdb content (optional). If not specified design an unconditional monomer
                    :param contig: contigs (around what and how to design a protein). See https://github.com/RosettaCommons/RFdiffusion
                    :param symmetry: To design a symmetrical protein,
                      Available symmetries:
                      - Cyclic symmetry (C_n) # call as c5
                      - Dihedral symmetry (D_n) # call as d5
                      - Tetrahedral symmetry # call as tetrahedral
                      - Octahedral symmetry # call as octahedral
                      - Icosahedral symmetry # call as icosahedral
                    Default None
                    :param timesteps: default 50 ( Desired iterations to generate structure. )
                    :param hotspots: A30, A33, A34
                    The model optionally readily learns that it should be making an interface which involving these hotspot residues. Input is ChainResidueNumber: A100 for residue 100 on chain A.
                    :return:
                    """
    allowed_symmetry_enum = ['cyclic', 'dihedral', 'tetrahedral', 'octahedral', 'icosahedral']
    if symmetry and symmetry not in ['tetrahedral', 'octahedral', 'icosahedral']:
        raise Exception(f'Symmetry must be one of {str(allowed_symmetry_enum)}')

    try:
        inference_path = os.path.join(rfdiffusion_dir, 'scripts', 'run_inference')
        program = ['python3.9', inference_path]

        program.append(f'inference.output_prefix={output_files_dir}')

        if pdb_content:
            program.append(f'inference.input_pdb={input_pdbs_dir}')

        program.append(f'contigmap.contigs=[{contig}]')

        if symmetry:
            program.append('--config-name')
            program.append('symmetry')
            program.append(f'inference.symmetry={symmetry}')

        if timesteps:
            program.append(f'diffuser.partial_T={timesteps}')

        if hotspots:
            program.append(f'ppi.hotspot_res=[{hotspots}]')

        res = subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        print(res.stdout.decode('utf-8'))
        print(res.stderr.decode('utf-8'))

        result = []
        files = glob.glob(os.path.join(output_files_dir, '*'))
        for f in files:
            with open(f, 'r') as f:
                result.append(f.read())
        return result
    except Exception:
        print(traceback.format_exc())
    finally:
        cleanup_folders()

    return []
