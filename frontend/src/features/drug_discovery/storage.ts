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
  getTargetData,
  getLigandData,
  predictFolding,
  predictBindingPocket
  // other API functions
} from 'src/features/drug_discovery/api';

import {
  Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
  Body_upload_target_api_v1_drug_discovery_upload_target_post,
  TargetMetaData,
  LigandMetaData, ExperimentMetadataResponse
} from 'src/api/client';

export interface TargetData {
  targetId: string;
  sequence: string;
  pdbContents: string | null;
}


export const useDrugDiscoveryStore = defineStore('drugDiscovery', {

  state: () => ({
    experiments: [] as ExperimentMetadataResponse[],
    currentExperiment: null as ExperimentMetadataResponse | null,
    targets: [] as TargetMetaData[],
    currentTarget: null,
    targetData: null as TargetData | null,
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
        this.experiments = this.experiments.filter(exp => exp.id !== experimentId);
      } catch (error) {
        console.error('Error deleting experiment:', error);
      }
    },
    async uploadTargetToExperiment(experimentId: string, targetFile: File) {
      try {
        const targetData: Body_upload_target_api_v1_drug_discovery_upload_target_post = {
          experimentId: experimentId,
          fastaFile: targetFile
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
          experimentId: experimentId,
          targetId: targetId,
          sdfFile: sdfFile
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
        this.targets = response.map(target => ({
          targetId: target.target_id,
          targetName: target.target_name
        }));
        // Optionally set the currentExperiment to the selected one
      } catch (error) {
        console.error('Error fetching targets:', error);
      }
    },
    async fetchLigandsForTarget(experimentId: string, targetId: string) {
      try {
        const response = await getLigandsList(experimentId, targetId);
        return response.map(ligand => ({
          ligandId: ligand.ligand_id,
          ligandName: ligand.ligand_name
        }));
        // Optionally set the currentTarget to the selected one
      } catch (error) {
        console.error('Error fetching ligands:', error);
      }
    },
    async fetchTargetData(experimentId: string, targetId: string) {
      try {
        const response = await getTargetData(experimentId, targetId);
        return {
          targetId: targetId,
          sequence: response.protein_sequence,
          pdbContents: response.protein_pdb || null,
        };
        // Optionally set the currentTarget to the selected one
      } catch (error) {
        console.error('Error fetching ligands:', error);
      }
    },
    async predictFoldingForTarget(experimentId: string, targetId: string) {
      try {
        const response = await predictFolding(experimentId, targetId);
        return response.pdb_content;
        // Optionally set the currentTarget to the selected one
      } catch (error) {
        console.error('Error fetching ligands:', error);
      }
    },
    async predictPocketForTarget(experimentId: string, targetId: string) {
      try {
        const response = await predictBindingPocket(experimentId, targetId);
        return response.pocket_ids;
        // Optionally set the currentTarget to the selected one
      } catch (error) {
        console.error('Error fetching ligands:', error);
      }
    },
    async fetchLigandData(experimentId: string, targetId: string, ligandId: string) {
      try {
        const response = await getLigandData(experimentId, targetId, ligandId);
        return {
          ligandId: ligandId,
          ligandName: response.ligand_name,
          ligandSdf: response.ligand_sdf,
          ligandSmiles: response.ligand_smiles
        };
        // Optionally set the currentTarget to the selected one
      } catch (error) {
        console.error('Error fetching ligands:', error);
      }
    },
  }
});