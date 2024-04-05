import os.path

__all__ = ['Settings', ]

from nolabs.infrastructure import environment
from nolabs.infrastructure.environment import is_dev

import configparser
import os
import nolabs

class Settings:
    def __init__(self):
        self._config = configparser.ConfigParser()
        dir_path = os.path.dirname(os.path.realpath(__file__))

        settings_path = 'settings.ini' if is_dev() else 'settings.ini'

        self._config.read(os.path.join(dir_path, settings_path))

    @property
    def is_test(self) -> bool:
        return self._config.getboolean('modules', 'test')

    @property
    def solubility_host(self) -> str:
        host = self._config.get('microservices', 'solubility')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'solubility')
        return host

    @property
    def conformations_host(self) -> str:
        host = self._config.get('microservices', 'conformations')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'conformations')
        return host

    @property
    def localisation_host(self) -> str:
        host = self._config.get('microservices', 'localisation')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'localisation')
        return host

    @property
    def gene_ontology_host(self) -> str:
        host = self._config.get('microservices', 'gene_ontology')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'gene_ontology')
        return host

    @property
    def p2rank_host(self) -> str:
        host = self._config.get('microservices', 'p2rank')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'p2rank')
        return host

    @property
    def msa_light_host(self) -> str:
        host = self._config.get('microservices', 'msa_light')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'msa_light')
        return host

    @property
    def umol_host(self) -> str:
        host = self._config.get('microservices', 'umol')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'umol')
        return host

    @property
    def diffdock_host(self) -> str:
        return self._config.get('microservices', 'diffdock')

    @property
    def biobuddy_host(self) -> str:
        return self._config.get('microservices', 'biobuddy')

    @property
    def rcsb_pdb_query_host(self) -> str:
        return self._config.get('microservices', 'rcsb_pdb_query')

    @property
    def pubmed_query_host(self) -> str:
        return self._config.get('microservices', 'pubmed_query')

    @property
    def chembl_query_host(self) -> str:
        return self._config.get('microservices', 'chembl_query')

    @property
    def esmfold_light_host(self) -> str:
        host = self._config.get('microservices', 'esmfold_light')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'esmfold_light')
        return host

    @property
    def rosettafold_host(self) -> str:
        host = self._config.get('microservices', 'rosettafold')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'rosettafold')
        return host

    @property
    def esmfold_host(self) -> str:
        host = self._config.get('microservices', 'esmfold')
        if environment.is_compose():
            return host.replace('127.0.0.1', 'esmfold')
        return host

    @property
    def conformations_simulations_file_name(self) -> str:
        return self._config.get('conformations', 'file_name')

    @property
    def localisation_experiments_folder(self) -> str:
        import nolabs
        exp_folder = self._config.get('localisation', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def conformations_experiments_folder(self) -> str:
        import nolabs
        exp_folder = self._config.get('conformations', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def conformations_metadata_file_name(self) -> str:
        return self._config.get('conformations', 'metadata_file_name')

    @property
    def small_molecules_experiments_folder(self) -> str:
        exp_folder = self._config.get('small-molecules-design', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def small_molecules_metadata_file_name(self) -> str:
        return self._config.get('small-molecules-design', 'metadata_file_name')


    @property
    def protein_design_experiments_folder(self) -> str:
        exp_folder = self._config.get('protein-design', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def protein_design_host(self) -> str:
        return self._config.get('microservices', 'protein_design')

    @property
    def protein_design_metadata_file_name(self) -> str:
        return self._config.get('protein-design', 'metadata_file_name')

    @property
    def go_file_name(self) -> str:
        return self._config.get('go', 'file_name')

    @property
    def folding_experiments_folder(self) -> str:
        exp_folder = self._config.get('folding', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def folding_metadata_file_name(self) -> str:
        return self._config.get('folding', 'metadata_file_name')

    @property
    def go_experiments_folder(self) -> str:
        exp_folder = self._config.get('go', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def go_metadata_file_name(self) -> str:
        return self._config.get('go', 'metadata_file_name')

    @property
    def solubility_experiments_folder(self) -> str:
        exp_folder = self._config.get('solubility', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def solubility_metadata_file_name(self) -> str:
        return self._config.get('solubility', 'metadata_file_name')

    @property
    def drug_discovery_experiments_folder(self):
        exp_folder = self._config.get('drug-discovery', 'experiments_path')
        return os.path.join(os.path.dirname(nolabs.__file__), exp_folder)

    @property
    def drug_discovery_experiment_metadata_file_name(self) -> str:
        return self._config.get('drug-discovery', 'experiment_metadata_file_name')

    @property
    def drug_discovery_target_metadata_file_name(self) -> str:
        return self._config.get('drug-discovery', 'target_metadata_file_name')

    @property
    def drug_discovery_ligand_metadata_file_name(self) -> str:
        return self._config.get('drug-discovery', 'ligand_metadata_file_name')

    @property
    def drug_discovery_pocket_directory_name(self) -> str:
        return self._config.get('drug-discovery', 'pocket_directory_name')

    @property
    def drug_discovery_pocket_file_name(self) -> str:
        return self._config.get('drug-discovery', 'pocket_file_name')

    @property
    def folding_host(self) -> str:
        return self._config.get('microservices', 'esmfold_host')

    @property
    def drug_discovery_self_hosted_msa(self) -> str:
        return self._config.get('drug-discovery', 'self_hosted_msa')

    @property
    def drug_discovery_msa_remote_prediction_url(self) -> str:
        return self._config.get('drug-discovery', 'msa_server_url')
    @property
    def drug_discovery_msa_file_name(self) -> str:
        return self._config.get('drug-discovery', 'msa_file_name')
    @property
    def drug_discovery_running_jobs_dir_name(self) -> str:
        return self._config.get('drug-discovery', 'running_job_dir_name')

    @property
    def drug_discovery_running_jobs_config_file_name(self) -> str:
        return self._config.get('drug-discovery', 'running_jobs_config_name')

    @property
    def drug_discovery_running_jobs_metadata_file_name(self) -> str:
        return self._config.get('drug-discovery', 'jobs_metadata_file_name')

    @property
    def drug_discovery_running_jobs_pocket_file_name(self) -> str:
        return self._config.get('drug-discovery', 'job_pocket_file_name')

    @property
    def drug_discovery_running_jobs_progress_file_name(self) -> str:
        return self._config.get('drug-discovery', 'job_progress_file')

    @property
    def drug_discovery_docking_result_directory_name(self) -> str:
        return self._config.get('drug-discovery', 'docking_result_directory_name')

    @property
    def drug_discovery_docking_result_metadata_filename_name(self) -> str:
        return self._config.get('drug-discovery', 'docking_result_metadata_file_name')

    @property
    def drug_discovery_umol_docking_result_sdf_file_name(self) -> str:
        return self._config.get('drug-discovery', 'umol_docking_result_sdf_file_name')

    @property
    def drug_discovery_umol_docking_result_pdb_file_name(self) -> str:
        return self._config.get('drug-discovery', 'umol_docking_result_pdb_file_name')

    @property
    def drug_discovery_umol_docking_result_plddt_file_name(self) -> str:
        return self._config.get('drug-discovery', 'umol_docking_result_plddt_file_name')

    @property
    def drug_discovery_diffdock_docking_result_pdb_file_name(self) -> str:
        return self._config.get('drug-discovery', 'diffdock_docking_result_pdb_file_name')

    @property
    def drug_discovery_diffdock_docking_results_metadata_file_name(self) -> str:
        return self._config.get('drug-discovery', 'diffdock_docking_results_metadata_file_name')

    @property
    def drug_discovery_diffdock_params_file_name(self) -> str:
        return self._config.get('drug-discovery', 'diffdock_params_file_name')


    @property
    def biobuddy_conversation_file_name(self) -> str:
        return self._config.get('biobuddy', 'conversation_file_name')

    @property
    def reinvent_host(self) -> str:
        return self._config.get('microservices', 'reinvent')