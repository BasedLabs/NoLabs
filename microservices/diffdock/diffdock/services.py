import tempfile
import subprocess
import os
import csv
import glob
import re
import shutil
from fastapi import HTTPException
from diffdock.api_models import RunDiffDockPredictionRequest, RunDiffDockPredictionResponse, SDFResult

__all__ = ['run_docking']

os.environ["OMP_NUM_THREADS"] = "4"  # Example: Use 4 threads. Adjust based on your CPU.
os.environ["MKL_NUM_THREADS"] = "4"

def clear_results_directory(directory):
    # Check if the directory exists
    if os.path.exists(directory):
        # Remove all files and subdirectories in the directory
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            if os.path.isfile(item_path) or os.path.islink(item_path):
                os.unlink(item_path)  # Remove files and links
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)  # Recursively remove subdirectories

def run_docking(request: RunDiffDockPredictionRequest) -> RunDiffDockPredictionResponse:

    results_dir = "/app/DiffDock/results/user_predictions_small_new"
    
    clear_results_directory(results_dir)

    # Create temporary files for the protein and ligand
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", mode="w+") as protein_file, \
            tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", mode="w+") as ligand_file, \
            tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode="w+", newline='') as csv_file:
        
        protein_file.write(request.pdb_contents)
        ligand_file.write(request.sdf_contents)
        protein_file.flush()
        ligand_file.flush()

        # Write the protein and ligand paths to a CSV file
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["protein_path", "ligand"])
        csv_writer.writerow([protein_file.name, ligand_file.name])
        csv_file.flush()

        subprocess_env = os.environ.copy()
        subprocess_env["OMP_NUM_THREADS"] = "4"  # Adjust based on your CPU
        subprocess_env["MKL_NUM_THREADS"] = "4"


        # Prepare fasta file
        print("preparing fasta file")
        subprocess.run(["python", "datasets/esm_embedding_preparation.py", "--protein_ligand_csv", csv_file.name, "--out_file", "data/prepared_for_esm.fasta"], check=True, cwd="/app/DiffDock", env=subprocess_env)

        # Generate ESM embeddings
        print("Generating embeddings")
        os.environ['HOME'] = '/app/DiffDock/esm/model_weights'
        os.environ['PYTHONPATH'] = os.environ.get('PYTHONPATH', '') + ':/app/DiffDock/esm'
        subprocess.run(["python", "/app/DiffDock/esm/scripts/extract.py", "esm2_t33_650M_UR50D", "data/prepared_for_esm.fasta", "data/esm2_output", "--repr_layers", "33", "--include", "per_tok", "--truncation_seq_length", "30000"], check=True, cwd="/app/DiffDock", env=subprocess_env)

        # Run DiffDock inference
        print("Running inference")
        result = subprocess.run(["python", "-m", "inference", "--protein_ligand_csv", csv_file.name, "--out_dir", "results/user_predictions_small_new", "--inference_steps", str(request.inference_steps), "--samples_per_complex", str(request.samples_per_complex), "--batch_size", str(request.batch_size)], capture_output=True, text=True, cwd="/app/DiffDock", env=subprocess_env)

        if result.returncode != 0:
            raise Exception(result.stderr)

        # Process the results
        results_dir = "/app/DiffDock/results/user_predictions_small_new"
        sdf_results = []
        
        # Read the original PDB file contents
        with open(protein_file.name, 'r') as file:
            pdb_contents = file.read()

        # List all .sdf files in the results directory
        for result_dir in glob.glob(os.path.join(results_dir, "index*")):
            if os.path.isdir(result_dir):
                for sdf_file in glob.glob(os.path.join(result_dir, "*.sdf")):
                    with open(sdf_file, 'r') as file:
                        sdf_content = file.read()

                    # Extract confidence from the file name
                    confidence_match = re.search("confidence([\-\.\d]+)\.sdf", sdf_file)
                    if confidence_match:
                        confidence = float(confidence_match.group(1))
                    else:
                        confidence = None  # Handle the case where confidence is not in the file name

                    # Run gnina to get scored affinity
                    scored_stdout = subprocess.check_output(["/app/DiffDock/gnina", "--score_only", "-r", protein_file.name, "-l", sdf_file])
                    scored_affinity_match = re.search("Affinity:\s*([\-\.\d+]+)", scored_stdout.decode('utf-8'))
                    scored_affinity = float(scored_affinity_match.group(1)) if scored_affinity_match else None

                    # Run gnina to get minimized affinity
                    minimized_stdout = subprocess.check_output(["/app/DiffDock/gnina", "--local_only", "--minimize", "-r", protein_file.name, "-l", sdf_file, "--autobox_ligand", sdf_file, "--autobox_add", "2"])
                    minimized_affinity_match = re.search("Affinity:\s*([\-\.\d+]+)", minimized_stdout.decode('utf-8'))
                    minimized_affinity = float(minimized_affinity_match.group(1)) if minimized_affinity_match else None

                    # Append the result for this SDF file
                    sdf_results.append(SDFResult(
                        sdf_content=sdf_content,
                        confidence=confidence,
                        scored_affinity=scored_affinity,
                        minimized_affinity=minimized_affinity
                    ))

        # Construct the response
        response = RunDiffDockPredictionResponse(
            success=True,
            message="Docking completed successfully",
            pdb_contents=pdb_contents,
            sdf_results=sdf_results
        )

    return response