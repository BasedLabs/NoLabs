import glob
import os
import shutil
import subprocess

from protein_design.api_models import (RunRfdiffusionRequest,
                                       RunRfdiffusionResponse)
from protein_design.loggers import logger
from protein_design.shared import JobStatusEnum, memory_manager


def run_rfdiffusion(request: RunRfdiffusionRequest) -> RunRfdiffusionResponse:
    memory_manager.change_status(str(request.job_id), JobStatusEnum.running)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    rfdiffusion_dir = os.path.join(current_dir, "RFdiffusion")
    input_pdbs_dir = "input_pdbs"
    output_files_dir = "output_pdbs"
    tmp_pdb_file = "tmp.pdb"

    if not os.path.exists(input_pdbs_dir):
        os.mkdir(input_pdbs_dir)

    if not os.path.exists(output_files_dir):
        os.mkdir(output_files_dir)

    try:
        inference_path = os.path.join(rfdiffusion_dir, "scripts", "run_inference.py")
        program = [
            "python3.9",
            inference_path,
            f"inference.model_directory_path={rfdiffusion_dir}/models",
            f"inference.num_designs={request.number_of_designs}",
            f"inference.output_prefix={output_files_dir}/result",
        ]

        if request.pdb_content:
            with open(tmp_pdb_file, "w") as f:
                f.write(request.pdb_content)
            program.append(f"inference.input_pdb={tmp_pdb_file}")

        program.append(f"contigmap.contigs=[{request.contig}]")

        program.append("denoiser.noise_scale_ca=0")
        program.append("denoiser.noise_scale_frame=0")
        program.append(f"diffuser.T={request.timesteps}")

        if request.hotspots:
            program.append(f"ppi.hotspot_res=[{request.hotspots}]")

        res = subprocess.run(
            program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False
        )

        stderr = res.stderr.decode("utf-8")

        pdbs = []
        files = glob.glob(os.path.join(output_files_dir, "*.pdb"))
        for file_name in files:
            with open(file_name, "r") as f:
                pdbs.append(f.read())
        logger.exception(stderr)
        return RunRfdiffusionResponse(
            pdbs_content=pdbs,
            errors=[f"Unknown error, try to fix contig or hotspots input, {stderr}"],
        )
    except Exception as e:
        logger.exception(str(e))
        return RunRfdiffusionResponse(
            pdbs_content=[], errors=["Internal exception occured"]
        )
    finally:
        shutil.rmtree(input_pdbs_dir)
        shutil.rmtree(output_files_dir)
        memory_manager.change_status(str(request.job_id), JobStatusEnum.idle)
