import dataclasses
import json
import os.path
import shutil
from typing import List, Tuple

import numpy as np

from nolabs.domain.experiment import ExperimentId
from nolabs.infrastructure.settings import Settings
from nolabs.utils.pdb import PDBWriter, PDBReader

from nolabs.utils.sdf import SDFReader, SDFWriter
from nolabs.modules.drug_discovery.data_models.target import TargetId
from nolabs.modules.drug_discovery.data_models.ligand import LigandId
from nolabs.modules.drug_discovery.data_models.result import JobId, UmolResultMetaData, \
    UmolDockingResultData, JobMetaData, DiffDockDockingResultData, DiffDockResultMetaData, DiffDockLigandResultData, \
    DiffDockModelParams
from nolabs.modules.drug_discovery.services.ligand_file_management import LigandsFileManagement


class ResultsFileManagement:
    def __init__(self, settings: Settings, ligand_file_management: LigandsFileManagement):
        self._settings = settings
        self.sdf_reader = SDFReader()
        self.sdf_writer = SDFWriter()
        self.pdb_writer = PDBWriter()
        self.pdb_reader = PDBReader()
        self._ligand_file_management = ligand_file_management

    def ensure_results_folder_exists(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId):
        if not os.path.isdir(self.results_folder(experiment_id, target_id, ligand_id)):
            os.mkdir(self.results_folder(experiment_id, target_id, ligand_id))

    def ensure_result_folder_exists(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                                    job_id: JobId):
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        if not os.path.isdir(result_folder):
            os.mkdir(result_folder)

    def results_folder(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId) -> str:
        parent_ligand_folder = self._ligand_file_management.target_ligand_folder(experiment_id, target_id, ligand_id)
        return os.path.join(parent_ligand_folder, 'results')

    def result_folder(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                      job_id: JobId) -> str:
        return os.path.join(self.results_folder(experiment_id, target_id, ligand_id), job_id.value)

    def create_result_folder(
            self, experiment_id: ExperimentId,
            target_id: TargetId,
            ligand_id: LigandId,
            job_id: JobId
    ):
        self.ensure_results_folder_exists(experiment_id, target_id, ligand_id)
        results_dir = self.results_folder(experiment_id, target_id, ligand_id)
        if not os.path.exists(os.path.join(results_dir, job_id.value)):
            os.mkdir(os.path.join(results_dir, job_id.value))
            self.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "target_id", target_id.value)
            self.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "ligand_id", ligand_id.value)
            self.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "job_id", job_id.value)

    def update_job_folding_method(self, experiment_id: ExperimentId,
                                  target_id: TargetId,
                                  ligand_id: LigandId,
                                  job_id: JobId, folding_method: str):
        self.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "folding_method", folding_method)

    def update_job_docking_method(self, experiment_id: ExperimentId,
                                  target_id: TargetId,
                                  ligand_id: LigandId,
                                  job_id: JobId, docking_method: str):
        self.update_job_metadata(experiment_id, target_id, ligand_id, job_id, "docking_method", docking_method)

    def store_result_input_pocketIds(self,
                                     experiment_id: ExperimentId,
                                     target_id: TargetId,
                                     ligand_id: LigandId,
                                     job_id: JobId,
                                     pocket_ids: List[int]):
        self.ensure_results_folder_exists(experiment_id, target_id, ligand_id)
        result_dir = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        pocket_file = os.path.join(result_dir, self._settings.drug_discovery_running_jobs_pocket_file_name)
        pocket_arr = np.asarray(pocket_ids)
        np.save(pocket_file, pocket_arr)

    def get_result_input_pocketIds(self,
                                   experiment_id: ExperimentId,
                                   target_id: TargetId,
                                   ligand_id: LigandId,
                                   job_id: JobId) -> List[int]:
        self.ensure_results_folder_exists(experiment_id, target_id, ligand_id)
        result_dir = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        pocket_file = os.path.join(result_dir, self._settings.drug_discovery_running_jobs_pocket_file_name)
        return np.load(pocket_file)

    def check_binding_pocket_exist(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                                   job_id: JobId) -> bool:
        result_dir = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        pocket_file = os.path.join(result_dir, self._settings.drug_discovery_running_jobs_pocket_file_name)
        return os.path.exists(pocket_file)

    def store_umol_result_data(self, experiment_id: ExperimentId,
                               target_id: TargetId,
                               ligand_id: LigandId,
                               job_id: JobId,
                               result_data: UmolDockingResultData) -> UmolResultMetaData:

        self.create_result_folder(experiment_id, target_id, ligand_id, job_id)
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)

        sdf_file_name = self._settings.drug_discovery_umol_docking_result_sdf_file_name
        sdf_file_path = os.path.join(result_folder, sdf_file_name)
        sdf_contents = result_data.predicted_sdf
        self.sdf_writer.write_sdf(sdf_contents, sdf_file_path)

        pdb_file_name = self._settings.drug_discovery_umol_docking_result_pdb_file_name
        pdb_file_path = os.path.join(result_folder, pdb_file_name)
        pdb_contents = result_data.predicted_pdb
        self.pdb_writer.write_pdb(pdb_contents, pdb_file_path)

        plddt_file_name = self._settings.drug_discovery_umol_docking_result_plddt_file_name
        plddt_file_path = os.path.join(result_folder, plddt_file_name)
        plddt_list = result_data.plddt_array
        np.save(plddt_file_path, np.asarray(plddt_list))

        return UmolResultMetaData(job_id=job_id.value, target_id=target_id.value, ligand_id=ligand_id.value)

    def store_diffdock_result_data(self, experiment_id: ExperimentId,
                                   target_id: TargetId,
                                   ligand_id: LigandId,
                                   job_id: JobId,
                                   result_data: DiffDockDockingResultData) -> DiffDockResultMetaData:

        self.create_result_folder(experiment_id, target_id, ligand_id, job_id)
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)

        sdf_file_path = os.path.join(result_folder, result_data.predicted_sdf_file_name)
        self.sdf_writer.write_sdf(result_data.predicted_sdf_contents, sdf_file_path)

        pdb_file_name = self._settings.drug_discovery_diffdock_docking_result_pdb_file_name
        pdb_file_path = os.path.join(result_folder, pdb_file_name)
        pdb_contents = result_data.predicted_pdb
        if not os.path.exists(pdb_file_path):
            self.pdb_writer.write_pdb(pdb_contents, pdb_file_path)

        # Check if the JSON file already exists and has content
        json_file_path = os.path.join(result_folder,
                                      self._settings.drug_discovery_diffdock_docking_results_metadata_file_name)
        if os.path.exists(json_file_path) and os.path.getsize(json_file_path) > 0:
            # Load the existing content
            with open(json_file_path, 'r') as json_file:
                docking_results = json.load(json_file)
        else:
            docking_results = {}

        key = result_data.predicted_sdf_file_name  # Assume each result has a unique file name
        docking_results[key] = {
            "confidence": result_data.confidence,
            "scored_affinity": result_data.scored_affinity,
            "minimized_affinity": result_data.minimized_affinity
        }

        # Write the updated dictionary back to the JSON file
        with open(json_file_path, 'w') as json_file:
            json.dump(docking_results, json_file, indent=4)

        return DiffDockResultMetaData(job_id=job_id.value,
                                      target_id=target_id.value,
                                      ligand_id=ligand_id.value,
                                      scored_affinity=result_data.scored_affinity,
                                      minimized_affinity=result_data.minimized_affinity,
                                      confidence=result_data.confidence,
                                      ligand_file_name=result_data.predicted_sdf_file_name)

    def get_diffdock_ligand_data(self, experiment_id: ExperimentId, target_id: TargetId,
                                 ligand_id: LigandId,
                                 job_id: JobId,
                                 sdf_file_name: str) -> DiffDockLigandResultData:

        self.create_result_folder(experiment_id, target_id, ligand_id, job_id)
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        sdf_file_path = os.path.join(result_folder, sdf_file_name)
        sdf_contents = self.sdf_reader.read_sdf(sdf_file_path)

        return DiffDockLigandResultData(predicted_sdf_contents=sdf_contents)

    def delete_result(self, experiment_id: ExperimentId,
                      target_id: TargetId,
                      ligand_id: LigandId,
                      job_id: JobId) -> JobId:
        self.results_folder(experiment_id, target_id, ligand_id)
        self.result_folder(experiment_id, target_id, ligand_id, job_id)
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        shutil.rmtree(result_folder)

        return job_id

    def update_job_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                            job_id: JobId, key: str, value: str):
        metadata_file = os.path.join(self.result_folder(experiment_id, target_id, ligand_id, job_id),
                                     self._settings.drug_discovery_docking_result_metadata_filename_name)
        if os.path.exists(metadata_file):
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
                metadata[key] = value
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f)
        else:
            metadata = {key: value}
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f)

    def get_result_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                            job_id: JobId) -> UmolResultMetaData:
        metadata_file = os.path.join(self.result_folder(experiment_id, target_id, ligand_id, job_id),
                                     self._settings.drug_discovery_ligand_metadata_file_name)
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        result_metadata = UmolResultMetaData(job_id=metadata["job_id"],
                                             target_id=metadata["target_id"],
                                             ligand_id=metadata["ligand_id"])
        return result_metadata

    def get_job_metadata(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                         job_id: JobId) -> JobMetaData:
        metadata_file = os.path.join(self.result_folder(experiment_id, target_id, ligand_id, job_id),
                                     self._settings.drug_discovery_ligand_metadata_file_name)
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        result_metadata = JobMetaData(job_id=metadata["job_id"],
                                      target_id=metadata["target_id"],
                                      ligand_id=metadata["ligand_id"],
                                      folding_method=metadata["folding_method"],
                                      docking_method=metadata["docking_method"])
        return result_metadata

    def get_diffdock_params(self, experiment_id: ExperimentId,
                            target_id: TargetId,
                            ligand_id: LigandId,
                            job_id: JobId) -> DiffDockModelParams:
        params_file = os.path.join(self.result_folder(experiment_id, target_id, ligand_id, job_id),
                                   self._settings.drug_discovery_diffdock_params_file_name)

        if not os.path.exists(params_file):
            default_params = {
                "inference_steps": 20,
                "samples_per_complex": 1,
                "batch_size": 10,
                "actual_steps": 18
            }
            with open(params_file, 'w') as f:
                json.dump(default_params, f)

        with open(params_file, 'r') as f:
            params = json.load(f)

        return DiffDockModelParams(**params)

    def update_diffdock_params(self, experiment_id: ExperimentId, target_id: TargetId, ligand_id: LigandId,
                               job_id: JobId, new_params: DiffDockModelParams):
        params_file = os.path.join(self.result_folder(experiment_id, target_id, ligand_id, job_id),
                                   self._settings.drug_discovery_diffdock_params_file_name)

        with open(params_file, 'w') as f:
            # Convert the DiffDockModelParams instance to a dictionary and save it as JSON
            json.dump(dataclasses.asdict(new_params), f)

    def get_jobs_list(self, experiment_id: ExperimentId,
                      target_id: TargetId, ligand_id: LigandId) -> List[JobMetaData]:
        results_folder = self.results_folder(experiment_id, target_id, ligand_id)
        jobs_list = []
        if not os.path.exists(results_folder):
            return jobs_list
        for t_id in os.listdir(results_folder):
            if os.path.isdir(os.path.join(results_folder, t_id)):
                job_id = JobId(t_id)
                jobs_list.append(self.get_job_metadata(experiment_id, target_id, ligand_id, job_id))

        return jobs_list

    def get_umol_docking_result_data(self,
                                     experiment_id: ExperimentId,
                                     target_id: TargetId, ligand_id: LigandId,
                                     job_id: JobId) -> UmolDockingResultData:
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)

        sdf_file_name = self._settings.drug_discovery_umol_docking_result_sdf_file_name
        sdf_file_path = os.path.join(result_folder, sdf_file_name)
        sdf_contents = self.sdf_reader.read_sdf(sdf_file_path)

        pdb_file_name = self._settings.drug_discovery_umol_docking_result_pdb_file_name
        pdb_file_path = os.path.join(result_folder, pdb_file_name)
        pdb_contents = self.pdb_reader.read_pdb(pdb_file_path)

        plddt_file_name = self._settings.drug_discovery_umol_docking_result_plddt_file_name
        plddt_file_path = os.path.join(result_folder, plddt_file_name)
        plddt_list = np.load(plddt_file_path).tolist()

        return UmolDockingResultData(predicted_pdb=pdb_contents, predicted_sdf=sdf_contents, plddt_array=plddt_list)

    def check_umol_result_data_available(self, experiment_id: ExperimentId,
                                         target_id: TargetId,
                                         ligand_id: LigandId,
                                         job_id: JobId) -> bool:
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        if not os.path.exists(result_folder):
            return False

        sdf_file_name = self._settings.drug_discovery_umol_docking_result_sdf_file_name
        sdf_file_path = os.path.join(result_folder, sdf_file_name)
        if not os.path.exists(sdf_file_path):
            return False

        pdb_file_name = self._settings.drug_discovery_umol_docking_result_pdb_file_name
        pdb_file_path = os.path.join(result_folder, pdb_file_name)
        if not os.path.exists(pdb_file_path):
            return False

        plddt_file_name = self._settings.drug_discovery_umol_docking_result_plddt_file_name
        plddt_file_path = os.path.join(result_folder, plddt_file_name)
        if not os.path.exists(plddt_file_path):
            return False

        return True

    def check_diffdock_result_data_available(self, experiment_id: ExperimentId,
                                             target_id: TargetId,
                                             ligand_id: LigandId,
                                             job_id: JobId) -> bool:
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)
        if not os.path.exists(result_folder):
            return False

        pdb_file_name = self._settings.drug_discovery_diffdock_docking_result_pdb_file_name
        pdb_file_path = os.path.join(result_folder, pdb_file_name)
        if not os.path.exists(pdb_file_path):
            return False

        docking_metadata_file = self._settings.drug_discovery_diffdock_docking_results_metadata_file_name
        docking_metadata_file_path = os.path.join(result_folder, docking_metadata_file)
        if not os.path.exists(docking_metadata_file_path):
            return False

        return True

    def get_diffdock_docking_result_data(self,
                                         experiment_id: ExperimentId,
                                         target_id: TargetId, ligand_id: LigandId,
                                         job_id: JobId) -> Tuple[str, List[DiffDockResultMetaData]]:
        result_folder = self.result_folder(experiment_id, target_id, ligand_id, job_id)

        # Path to the JSON metadata file
        json_file_path = os.path.join(result_folder,
                                      self._settings.drug_discovery_diffdock_docking_results_metadata_file_name)

        # Path to the PDB file
        pdb_file_name = self._settings.drug_discovery_diffdock_docking_result_pdb_file_name
        pdb_file_path = os.path.join(result_folder, pdb_file_name)

        if not os.path.exists(json_file_path):
            raise FileNotFoundError(f"No results metadata found for the specified IDs at {json_file_path}")

        if not os.path.exists(pdb_file_path):
            raise FileNotFoundError(f"No PDB file found for the specified IDs at {pdb_file_path}")

        # Read the PDB file contents
        with open(pdb_file_path, 'r') as pdb_file:
            predicted_pdb = pdb_file.read()

        # Read the JSON metadata file
        with open(json_file_path, 'r') as json_file:
            docking_results = json.load(json_file)

        ligand_metadata_list = []
        for ligand_file_name, data in docking_results.items():
            ligand_metadata = DiffDockResultMetaData(
                job_id=job_id.value,
                target_id=target_id.value,
                ligand_id=ligand_id.value,
                ligand_file_name=ligand_file_name,
                minimized_affinity=data["minimized_affinity"],
                scored_affinity=data["scored_affinity"],
                confidence=data.get("confidence")  # Use .get() in case "confidence" key might be missing
            )
            ligand_metadata_list.append(ligand_metadata)

        return predicted_pdb, ligand_metadata_list
