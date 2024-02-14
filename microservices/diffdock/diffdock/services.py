import tempfile
import subprocess
from fastapi import HTTPException
from diffdock.api_models import RunDiffDockPredictionRequest, RunDiffDockPredictionResponse

__all__ = ['run_docking']


def run_docking(request: RunDiffDockPredictionRequest) -> RunDiffDockPredictionResponse:
    try:
        # Create temporary files for the protein and ligand
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", mode="w+") as protein_file, \
             tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", mode="w+") as ligand_file:
            protein_file.write(request.pdb_contents)
            ligand_file.write(request.sdf_contents)
            protein_file.flush()
            ligand_file.flush()

            # Construct the command
            cmd = [
                "python", "-m", "inference",
                "--protein_path", protein_file.name,
                "--ligand_path", ligand_file.name,
                "--inference_steps", str(request.inference_steps),
                "--samples_per_complex", str(request.samples_per_complex),
                "--batch_size", str(request.batch_size),
                "--actual_steps", str(request.actual_steps)
            ]
            if request.no_final_step_noise:
                cmd.append("--no_final_step_noise")

            # Execute the command
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise Exception(result.stderr)

            # Parse the output and construct the response
            # You'll need to adapt this part based on the actual output of your DiffDock command
            response = RunDiffDockPredictionResponse(success=True, message="Docking completed successfully")

            return response

    except Exception as e:
        # Log the exception or handle it as needed
        raise HTTPException(status_code=500, detail=str(e))

    finally:
        # Clean up temporary files
        protein_file.close()
        ligand_file.close()
        os.unlink(protein_file.name)
        os.unlink(ligand_file.name)