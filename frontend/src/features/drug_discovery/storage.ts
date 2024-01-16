import { defineStore } from 'pinia';
import {
  addExperiment,
  deleteExperiment,
  getExperiments,
  uploadTarget,
  deleteTarget,
  getTargetsList,
  uploadLigand,
  deleteLigand,
  getLigandsList,
  // other API functions
} from './api';

import {
  nolabs__api_models__drug_discovery__ExperimentMetadataResponse,
  Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
  Body_upload_target_api_v1_drug_discovery_upload_target_post,
  TargetMetaData,
  LigandMetaData
} from 'api/client';

export const useDrugDiscoveryStore = defineStore('drugDiscovery', {
  state: () => ({
    experiments: [] as nolabs__api_models__drug_discovery__ExperimentMetadataResponse[],
    currentExperiment: null as nolabs__api_models__drug_discovery__ExperimentMetadataResponse | null,
    targets: [] as TargetMetaData[],
    currentTarget: null,
    ligands: [] as LigandMetaData[],
    currentLigand: null,
  }),
  actions: {
    async fetchExperiments() {
      try {
        const response = await getExperiments();
        this.experiments = response;
      } catch (error) {
        console.error('Error fetching experiments:', error);
      }
    },
    async createExperiment() {
      try {
        const response = await addExperiment();
        this.experiments.push(response);
      } catch (error) {
        console.error('Error creating experiment:', error);
      }
    },
    async removeExperiment(experimentId: string) {
      try {
        await deleteExperiment(experimentId);
        this.experiments = this.experiments.filter(exp => exp.experiment_id !== experimentId);
      } catch (error) {
        console.error('Error deleting experiment:', error);
      }
    },
    async uploadTargetToExperiment(experimentId: string, targetFile: File) {
      try {
        const targetData: Body_upload_target_api_v1_drug_discovery_upload_target_post = {
          experiment_id: experimentId,
          fasta: targetFile
        };

        const response = await uploadTarget(targetData);
      } catch (error) {
        console.error('Error uploading target:', error);
      }
    },
    async deleteTargetFromExperiment(experimentId: string, targetId: string) {
      try {
        await deleteTarget(experimentId, targetId);
        // Update the targets list for the specific experiment
      } catch (error) {
        console.error('Error deleting target:', error);
      }
    },
    async uploadLigandToTarget(experimentId: string, targetId: string, sdfFile: File) {
      try {
        const ligandData: Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post = {
          experiment_id: experimentId,
          target_id: targetId,
          sdf_file: sdfFile
        };
        const response = await uploadLigand(ligandData);
        // Updatethe ligands list for the specific target
      } catch (error) {
        console.error('Error uploading ligand:', error);
      }
    },
    async deleteLigandFromTarget(experimentId: string, targetId: string, ligandId: string) {
      try {
        await deleteLigand(experimentId, targetId, ligandId);
        // Update the ligands list for the specific target
      } catch (error) {
        console.error('Error deleting ligand:', error);
      }
    },
    async fetchTargetsForExperiment(experimentId: string) {
      try {
        const response = await getTargetsList(experimentId);
        this.targets = response;
        // Optionally set the currentExperiment to the selected one
      } catch (error) {
        console.error('Error fetching targets:', error);
      }
    },
    async fetchLigandsForTarget(experimentId: string, targetId: string) {
      try {
        const response = await getLigandsList(experimentId, targetId);
        this.ligands = response;
        // Optionally set the currentTarget to the selected one
      } catch (error) {
        console.error('Error fetching ligands:', error);
      }
    },
  },
});