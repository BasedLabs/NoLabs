import requests
from msa_light.api_models import RunMsaPredictionRequest, RunMsaPredictionResponse
from msa_light.loggers import Log
from msa_light.job_state_manager import job_state_manager

__all__ = ['predict_msa_service']


def predict_msa_service(parameters: RunMsaPredictionRequest) -> RunMsaPredictionResponse:
    try:
        input_fasta = "/app/query.fasta"
        output_dir = "/app/output"
        tmp_dir = "/app/tmp"

        # Write the FASTA content to a file
        with open(input_fasta, "w") as f:
            f.write(parameters.fasta_contents)

        # Ensure the output and temporary directories exist
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(tmp_dir, exist_ok=True)

        # Run MMseqs2 search to generate the .a3m file
        mmseqs_cmd = [
            "mmseqs", "search", input_fasta, MMSEQS_DB, output_dir, tmp_dir,
            "--msa-format-output", "a3m"
        ]

        try:
            subprocess.run(mmseqs_cmd, check=True)
        except subprocess.CalledProcessError as e:
            # Clean up before returning the error
            cleanup_files([input_fasta], [output_dir, tmp_dir])
            raise HTTPException(status_code=500, detail=f"MMseqs2 failed with error: {e}")

        # Find and read the .a3m file (adjust the path based on the MMseqs2 output structure)
        a3m_file = os.path.join(output_dir, "query.a3m")

        if not os.path.exists(a3m_file):
            # Clean up before returning the error
            cleanup_files([input_fasta], [output_dir, tmp_dir])
            raise HTTPException(status_code=500, detail="Failed to generate .a3m file")

        # Read and return the contents of the .a3m file
        with open(a3m_file, "r") as result_file:
            a3m_content = result_file.read()

        return RunMsaPredictionResponse(msa_contents=a3m_content)
    except Exception as e:
        Log.exception()
        return RunMsaPredictionResponse(msa_contents=None)
