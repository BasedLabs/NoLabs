import subprocess
import sys
import uuid
from typing import Tuple

import log
import toml
from leaf import DirectoryObject, FileObject
from settings import settings


class Reinvent:
    def __init__(self):
        self.fs = DirectoryObject(settings.CONFIGS_DIRECTORY)
        self.logger = log.logger

    async def prepare_pdbqt(self, pdb: bytes) -> Tuple[bytes, str]:
        self.logger.info("Preparing PDBQT")
        tmp = self.fs.add_directory("tmp")
        pdb_file = tmp.add_file(str(uuid.uuid4()) + ".pdb")
        pdb_file.write_bytes(pdb)
        pdbqt = self.fs.add_directory("tmp").add_file(str(uuid.uuid4()) + ".pdbqt")
        prepare_pdbqt = self.fs.files.first_or_default(
            lambda o: o.name == "prepare_pdbqt.sh"
        )
        subprocess.run([prepare_pdbqt.full_path, pdb_file.full_path, pdbqt.full_path])
        self.logger.info("Done preparing PDBQT")
        return (pdbqt.read_bytes(), pdbqt.name)

    def run_reinforcement_learning(self, config_id: str):
        self.logger.info("Running reinforcement learning")
        directory = self.fs.directories.first_or_default(lambda x: x.name == config_id)
        rl_config = directory.files.first_or_default(lambda o: o.name == "RL.toml")
        log_file = directory.files.first_or_default(lambda o: o.name == "output.log")
        error_file = directory.files.first_or_default(lambda o: o.name == "error.log")

        run_reinforcement_learning_shell = directory.first_or_default(
            lambda o: o.name == "start_reinforcement_learning.sh"
        )

        process = subprocess.Popen(
            [
                run_reinforcement_learning_shell.full_path,
                rl_config.full_path,
                error_file.full_path,
                log_file.full_path,
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
        directory = self.fs.directories.first_or_default(lambda x: x.name == config_id)

        scoring_config = directory.files.first_or_default(
            lambda o: o.name == "Scoring.toml"
        )

        log_file = directory.files.first_or_default(lambda o: o.name == "output.log")
        error_file = directory.files.first_or_default(lambda o: o.name == "error.log")

        chkpt = directory.files.first_or_default(lambda o: o.name == "rl_direct.chkpt")
        sampling_config: FileObject = directory.files.first_or_default(
            lambda o: o.name == "Sampling.toml"
        )
        sampling_output = directory.add_file("sampling_direct.csv")

        sampling_config_toml = sampling_config.read(toml.load, mode="r")
        sampling_config_toml["parameters"]["output_file"] = sampling_output.full_path
        sampling_config_toml["parameters"]["model_file"] = chkpt.full_path
        sampling_config_toml["parameters"][
            "num_smiles"
        ] = number_of_molecules_to_generate

        sampling_config.write_string(toml.dumps(sampling_config_toml))

        directory.files.first_or_default(
            lambda o: o.name == "sampling_direct.csv"
        ).write_string("")
        directory.files.first_or_default(
            lambda o: o.name == "scoring_input.smi"
        ).write_string("")
        directory.files.first_or_default(
            lambda o: o.name == "scoring_direct.csv"
        ).write_string("")

        shell = directory.files.first_or_default(
            lambda o: o.name == "start_sampling.sh"
        )

        process = subprocess.Popen(
            [
                shell.full_path,
                sampling_config.full_path,
                scoring_config.full_path,
                error_file.full_path,
                log_file.full_path,
            ],
            stdout=sys.stdout,
            stderr=sys.stderr,
        )

        process.wait()

        self.logger.info("Reinforcement sampling finished")
