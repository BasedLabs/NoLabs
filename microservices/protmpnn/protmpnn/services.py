import tempfile
import subprocess
import os
import shutil

from protmpnn.api_models import RunProtMPNNPredictionRequest, RunProtMPNNPredictionResponse

__all__ = ['run_protmpnn']

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

def run_protmpnn(request: RunProtMPNNPredictionRequest) -> RunProtMPNNPredictionResponse:
    results_dir = "/app/ProtMPNN/outputs"

    clear_results_directory(results_dir)

    # Create temporary files for the protein
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", mode="w+") as protein_file:
        protein_file.write(request.pdb_contents)
        protein_file.flush()

        # Depending on whether it's a monomer or homomer, we need to prepare the inputs differently

        if request.is_homomer:
            # For homomers, we can run protein_mpnn_run.py directly with --pdb_path and --pdb_path_chains

            # Prepare output directory
            output_dir = os.path.join(results_dir, "homomer_designs")
            os.makedirs(output_dir, exist_ok=True)

            # Build the command
            cmd = [
                "python", "/app/ProtMPNN/protein_mpnn_run.py",
                "--pdb_path", protein_file.name,
                "--pdb_path_chains", " ".join(request.chains_to_design) if request.chains_to_design else "",
                "--out_folder", output_dir,
                "--num_seq_per_target", str(request.num_seq_per_target),
                "--sampling_temp", str(request.sampling_temp),
                "--seed", str(request.seed),
                "--batch_size", str(request.batch_size)
            ]

            # Add fixed positions if provided
            if request.fixed_positions:
                # Prepare the fixed positions file
                fixed_positions_file = os.path.join(output_dir, "fixed_positions.txt")
                with open(fixed_positions_file, "w") as f:
                    for chain_id, positions in request.fixed_positions.items():
                        positions_str = ",".join(map(str, positions))
                        f.write(f"{chain_id}:{positions_str}\n")
                cmd.extend(["--fixed_positions_file", fixed_positions_file])

            # Run the command
            result = subprocess.run(cmd, capture_output=True, text=True, cwd="/app/ProtMPNN")
            if result.returncode != 0:
                raise Exception(result.stderr)

        else:
            # For monomers, we need to parse multiple chains and run with --jsonl_path

            # Prepare output directory
            output_dir = os.path.join(results_dir, "monomer_designs")
            os.makedirs(output_dir, exist_ok=True)

            # Prepare the input PDBs directory
            input_pdbs_dir = os.path.join(output_dir, "pdbs")
            os.makedirs(input_pdbs_dir, exist_ok=True)
            # Copy the protein file to the input PDBs directory
            input_pdb_path = os.path.join(input_pdbs_dir, os.path.basename(protein_file.name))
            shutil.copy(protein_file.name, input_pdb_path)

            # Prepare the parsed chains file
            path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")
            parse_cmd = [
                "python", "/app/ProtMPNN/helper_scripts/parse_multiple_chains.py",
                "--input_path", input_pdbs_dir,
                "--output_path", path_for_parsed_chains
            ]
            parse_result = subprocess.run(parse_cmd, capture_output=True, text=True, cwd="/app/ProtMPNN")
            if parse_result.returncode != 0:
                raise Exception(parse_result.stderr)

            # Build the command
            cmd = [
                "python", "/app/ProtMPNN/protein_mpnn_run.py",
                "--jsonl_path", path_for_parsed_chains,
                "--out_folder", output_dir,
                "--num_seq_per_target", str(request.num_seq_per_target),
                "--sampling_temp", str(request.sampling_temp),
                "--seed", str(request.seed),
                "--batch_size", str(request.batch_size)
            ]

            # Add fixed positions if provided
            if request.fixed_positions:
                # Prepare the fixed positions file
                fixed_positions_file = os.path.join(output_dir, "fixed_positions.txt")
                with open(fixed_positions_file, "w") as f:
                    for chain_id, positions in request.fixed_positions.items():
                        positions_str = ",".join(map(str, positions))
                        f.write(f"{chain_id}:{positions_str}\n")
                cmd.extend(["--fixed_positions_file", fixed_positions_file])

            # Run the command
            result = subprocess.run(cmd, capture_output=True, text=True, cwd="/app/ProtMPNN")
            if result.returncode != 0:
                raise Exception(result.stderr)

    # Now collect the outputs
    # The outputs are fasta files in the output directory
    sequences = []
    fasta_contents = []

    # Look for fasta files in output_dir
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".fa") or file.endswith(".fasta"):
                fasta_file_path = os.path.join(root, file)
                with open(fasta_file_path, 'r') as f:
                    fasta_content = f.read()
                    fasta_contents.append(fasta_content)
                    # Also extract sequences from fasta
                    # For simplicity, we can extract sequences directly
                    lines = fasta_content.strip().split('\n')
                    seq_lines = [line for line in lines if not line.startswith('>')]
                    sequence = ''.join(seq_lines)
                    sequences.append(sequence)

    response = RunProtMPNNPredictionResponse(
        success=True,
        message="Protein design completed successfully",
        sequences=sequences,
        fasta_contents=fasta_contents
    )

    return response
