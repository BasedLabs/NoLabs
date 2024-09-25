import subprocess
import sys
import uuid
from pathlib import Path
from typing import Tuple

import log
import toml
from settings import settings


class Reinvent:
    def __init__(self):
        self.fs = Path(settings.configs_directory)
        self.logger = log.logger

    async def prepare_pdbqt(self, pdb: bytes) -> Tuple[bytes, str]:
        self.logger.info("Preparing PDBQT")
        tmp = Path("/tmp")
        tmp.mkdir(exist_ok=True)
        pdb_file = tmp / (str(uuid.uuid4()) + ".pdb")
        pdb_file.touch()
        pdb_file.write_bytes(pdb)
        pdbqt = tmp / (str(uuid.uuid4()) + ".pdbqt")
        pdbqt.touch()
        prepare_pdbqt = self.fs / "base_configs" / "prepare_pdbqt.sh"
        subprocess.run([prepare_pdbqt, pdb_file, pdbqt])
        self.logger.info("Done preparing PDBQT")
        return (pdbqt.read_bytes(), str(pdbqt))

    def run_reinforcement_learning(self, config_id: str):
        self.logger.info("Running reinforcement learning")
        directory = self.fs / "jobs_configs" / config_id
        rl_config = directory / "RL.toml"
        log_file = directory / "output.log"
        error_file = directory / "error.log"

        run_reinforcement_learning_shell = directory / "start_reinforcement_learning.sh"

        process = subprocess.Popen(
            [
                run_reinforcement_learning_shell,
                rl_config,
                error_file,
                log_file,
            ],
            stdout=sys.stdout,
            stderr=sys.stderr,
        )

        process.wait()
        self.logger.info("Reinforcement learning finished")

    def run_sampling_and_scoring(
        self, config_id: str, number_of_molecules_to_generate: int
    ):
        self.logger.info("Running reinforcement sampling")
        directory = self.fs / "jobs_configs" / config_id
        scoring_config = directory / "Scoring.toml"
        log_file = directory / "output.log"
        error_file = directory / "error.log"
        directory / "rl_direct.chkpt"
        sampling_config: Path = directory / "Sampling.toml"
        sampling_output = directory / "sampling_direct.csv"
        sampling_output.touch()

        sampling_config_toml = toml.load(sampling_config.open("r"))
        sampling_config_toml["parameters"][
            "num_smiles"
        ] = number_of_molecules_to_generate
        sampling_config.write_text(toml.dumps(sampling_config_toml), encoding="utf-8")
        (directory / "scoring_direct.csv").touch()
        shell = directory / "start_sampling.sh"

        process = subprocess.Popen(
            [
                shell,
                sampling_config,
                scoring_config,
                error_file,
                log_file,
            ],
            stdout=sys.stdout,
            stderr=sys.stderr,
        )

        process.wait()

        self.logger.info("Reinforcement sampling finished")
