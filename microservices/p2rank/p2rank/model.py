import os
import subprocess
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from fastapi import UploadFile
from numpy import dtype, ndarray


class PocketPredictor:
    def __init__(self):
        pass

    # Method to get raw model outputs
    def _raw_inference(self, protein_file_path: str, save_dir: str) -> List[int]:
        protein_filename = os.path.split(protein_file_path)[1]
        ds = f"{save_dir}/protein_list.ds"
        with open(ds, "w") as out:
            out.write(f"{protein_filename}\n")
        p2rank_exec = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "p2rank_source",
            "p2rank_2.4.1",
            "prank",
        )
        cmd = [
            "bash",
            p2rank_exec,
            "predict",
            ds,
            "-o",
            f"{save_dir}/p2rank",
            "-threads",
            "1",
        ]

        # Run the command and wait for it to complete
        subprocess.run(cmd, check=True)

        p2rankFile = os.path.join(
            save_dir, "p2rank", f"{Path(protein_filename).stem}.pdb_predictions.csv"
        )
        pocket = pd.read_csv(p2rankFile, skipinitialspace=True)

        if pocket.empty:
            return []

        residue_ids = pocket["residue_ids"].str.split()
        all_ids = [
            int(item.split("_")[1]) for sublist in residue_ids for item in sublist
        ]

        pocket_ids = np.sort(np.array(all_ids))

        np.save(os.path.join(save_dir, "pocket.npy"), pocket_ids)

        pocket_ids = pocket_ids.tolist()

        return pocket_ids

    # Method to return raw outputs in the desired format
    def predict(self, pdb_contents: str) -> List[int]:
        temp_directory = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "temp"
        )

        if not os.path.exists(temp_directory):
            os.makedirs(temp_directory)

        temp_protein_pdb_path = os.path.join(temp_directory, "protein.pdb")

        with open(temp_protein_pdb_path, "w") as file:
            file.write(pdb_contents)

        return self._raw_inference(
            protein_file_path=temp_protein_pdb_path, save_dir=temp_directory
        )
