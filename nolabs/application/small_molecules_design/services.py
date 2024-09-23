import csv
import json
import shutil
import uuid
from distutils.dir_util import copy_tree
from pathlib import Path
from typing import List

import toml
from pydantic import BaseModel

from domain.models.small_molecules_design import SmallMoleculesDesignJob
from microservices.reinvent.service.api_models import PreparePdbqtRequest
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.settings import settings


class ReinventSmiles(BaseModel):
    smiles: str
    drug_likeness: float
    score: float
    stage: str


class ReinventSmilesRetriever:
    def __init__(self):
        self.reinvent_directory = settings.reinvent_directory

    def get_ligands(self, job_id: uuid.UUID) -> List[ReinventSmiles]:
        job_dir = self.reinvent_directory / str(job_id)

        rl_direct = job_dir / "rl_direct_1.csv"

        smiles = []

        if rl_direct.exists() and rl_direct.stat().st_size > 0:
            with rl_direct.open() as f:
                reader = csv.DictReader(f)
                smiles += list(reader)

        scoring_direct = job_dir / "scoring_direct.csv"

        if scoring_direct.exists() and scoring_direct.stat().st_size > 0:
            with scoring_direct.open() as f:
                reader = csv.DictReader(f)
                smiles += list(reader)

        for row in smiles:
            smiles = row["SMILES"]
            drug_likeness = row["QED"]
            score = row["Score"]
            stage = row["step"]

            smiles.append(
                ReinventSmiles(
                    smiles=smiles, drug_likeness=drug_likeness, score=score, stage=stage
                )
            )

        return smiles


class ReinventParametersSaver:
    def __init__(self):
        self.reinvent_directory = settings.reinvent_directory

    def _prepare_dockstream_config(
            self, job: SmallMoleculesDesignJob, pdbqt: Path, logfile: Path
    ) -> Path:
        # Configure dockstream
        dockstream_config = self.reinvent_directory / str(job.id) / "dockstream.json"

        dockstream_config_text = dockstream_config.read_text()
        dockstream_config_json = json.loads(dockstream_config_text)
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "parallelization"
        ][
            "number_cores"
        ] = 2  # type: ignore
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "receptor_pdbqt_path"
        ].append(pdbqt)
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "number_poses"
        ] = 2

        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "search_space"
        ]["--center_x"] = job.center_x
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "search_space"
        ]["--center_y"] = job.center_y
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "search_space"
        ]["--center_z"] = job.center_z
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "search_space"
        ]["--size_x"] = job.size_x
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "search_space"
        ]["--size_y"] = job.size_y
        dockstream_config_json["docking"]["docking_runs"][0]["parameters"][
            "search_space"
        ]["--size_z"] = job.size_z
        dockstream_config_json["docking"]["header"]["logging"]["logfile"] = logfile

        poses = self.reinvent_directory / str(job.id) / "poses.sdf"
        scores = self.reinvent_directory / str(job.id) / "scores.sdf"

        poses.touch()
        scores.touch()

        dockstream_config_json["docking"]["docking_runs"][0]["output"]["poses"][
            "poses_path"
        ] = poses
        dockstream_config_json["docking"]["docking_runs"][0]["output"]["scores"][
            "scores_path"
        ] = scores

        dockstream_config.write_text(json.dumps(dockstream_config_json, default=lambda x: str(x)))

        return dockstream_config

    def _prepare_scoring(self, job_dir: Path, dockstream_config: Path) -> Path:
        scoring_input = job_dir / "scoring_input.smi"
        scoring_input.touch()

        scoring_config = job_dir / "Scoring.toml"
        scoring_output = job_dir / "scoring_direct.csv"
        scoring_output.touch()

        scoring_config_toml = toml.loads(scoring_config.read_text())
        scoring_config_toml["output_csv"] = scoring_output
        scoring_config_toml["parameters"]["smiles_file"] = scoring_input
        scoring_config_toml["scoring"]["component"][0]["DockStream"]["endpoint"][0][
            "params"
        ]["configuration_path"] = str(dockstream_config)

        scoring_config.write_text(toml.dumps(scoring_config_toml))
        return scoring_config

    def _prepare_sampling(
            self, number_of_molecules_to_generate: int, job_dir: Path
    ) -> Path:
        chkpt = job_dir / "rl_direct.chkpt"
        sampling_config: Path = job_dir / "Sampling.toml"
        sampling_output = job_dir / "sampling_direct.csv"
        sampling_output.touch()

        sampling_config_toml = toml.loads(sampling_config.read_text())
        sampling_config_toml["parameters"]["output_file"] = sampling_output
        sampling_config_toml["parameters"]["model_file"] = chkpt
        sampling_config_toml["parameters"][
            "num_smiles"
        ] = number_of_molecules_to_generate

        sampling_config.write_text(toml.dumps(sampling_config_toml))
        return sampling_config

    def _prepare_rl(
            self,
            job_dir: Path,
            job: SmallMoleculesDesignJob,
            csv_result: Path,
            dockstream_config: Path,
    ) -> Path:
        # Configure learning
        reinforcement_learning_config = job_dir / "RL.toml"
        chkpt_file = job_dir / "rl_direct.chkpt"
        chkpt_file.touch()
        t = toml.loads(reinforcement_learning_config.read_text())
        t["parameters"]["batch_size"] = job.batch_size
        t["parameters"]["summary_csv_prefix"] = csv_result
        t["diversity_filter"]["minscore"] = job.minscore
        t["stage"][0]["chkpt_file"] = chkpt_file
        t["stage"][0]["max_steps"] = job.epochs
        t["stage"][0]["scoring"]["component"][0]["DockStream"]["endpoint"][0]["params"][
            "configuration_path"
        ] = dockstream_config

        reinforcement_learning_config.write_text(toml.dumps(t))
        return reinforcement_learning_config

    async def save_params(
            self, job: SmallMoleculesDesignJob, pdb: bytes | str
    ):
        job_dir = settings.reinvent_directory / str(job.id)
        if job_dir.exists():
            shutil.rmtree(job_dir)
        job_dir.mkdir(exist_ok=True, parents=True)
        copy_tree(str(Path("microservices") / "reinvent" / "reinvent_configs" / "base_configs"), str(job_dir))
        response = await celery.reinvent_prepare_target(
            task_id=uuid.uuid4(), payload=PreparePdbqtRequest(pdb=pdb)
        )
        csv_result = job_dir / "rl_direct"
        csv_result.touch()
        docking_log_file = job_dir / "docking.log"
        docking_log_file.touch()

        dockstream_config = self._prepare_dockstream_config(
            job=job, pdbqt=response.file_path, logfile=docking_log_file
        )
        self._prepare_rl(job_dir, job, csv_result, dockstream_config)
        self._prepare_sampling(2, job_dir)
        self._prepare_scoring(job_dir, dockstream_config)

        (job_dir / "output.log").touch()
        (job_dir / "error_log").touch()
