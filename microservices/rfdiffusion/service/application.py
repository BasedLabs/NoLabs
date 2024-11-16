import glob
import os
import subprocess
import logging

from settings import settings

from api_models import RunRfdiffusionRequest, RunRfdiffusionResponse


def design(request: RunRfdiffusionRequest) -> RunRfdiffusionResponse:
    logger = logging.getLogger(__name__)
    rfdiffusion_dir = settings.rfdiffusion_path
    input_pdbs_dir = 'input_pdbs'
    output_files_dir = 'output_pdbs'
    tmp_pdb_file = 'tmp.pdb'

    if not os.path.exists(input_pdbs_dir):
        os.mkdir(input_pdbs_dir)

    if not os.path.exists(output_files_dir):
        os.mkdir(output_files_dir)

    inference_path = os.path.join(rfdiffusion_dir, 'scripts', 'run_inference.py')
    program = ['python3.9', inference_path, f'inference.model_directory_path={rfdiffusion_dir}/models',
               f'inference.num_designs={request.number_of_designs}',
               f'inference.output_prefix={output_files_dir}/result']

    if request.pdb_content:
        with open(tmp_pdb_file, 'w') as f:
            f.write(request.pdb_content)
        program.append(f'inference.input_pdb={tmp_pdb_file}')

    if '[' in request.contig:
        request.contig = request.contig.lstrip('[').rstrip(']')

    program.append(f'contigmap.contigs=[{request.contig}]')

    program.append('denoiser.noise_scale_ca=0')
    program.append('denoiser.noise_scale_frame=0')
    program.append(f'diffuser.T={request.timesteps}')

    if request.hotspots:
        program.append(f'ppi.hotspot_res=[{request.hotspots}]')

    res = subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)

    stdout = res.stdout.decode('utf-8')
    stderr = res.stderr.decode('utf-8')

    logger.info(stdout)
    logger.error(stderr)

    pdbs = []
    files = glob.glob(os.path.join(output_files_dir, '*.pdb'))
    for file_name in files:
        with open(file_name, 'r') as f:
            pdbs.append(f.read())
    return RunRfdiffusionResponse(pdbs_content=pdbs, errors=[f'Unknown error, try to fix contig or hotspots input, {stderr}'])
