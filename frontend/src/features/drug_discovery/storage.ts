import {defineStore} from 'pinia';
import {
    createExperiment,
    deleteExperiment,
    deleteLigand,
    deleteTarget,
    getExperiments,
    getLigandData,
    getLigandsList,
    getTargetData,
    getTargetsList,
    predictBindingPocket,
    getTargetBindingPocket,
    predictFolding,
    uploadLigand,
    uploadTarget
} from 'src/features/drug_discovery/api';

import {
    Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
    Body_upload_target_api_v1_drug_discovery_upload_target_post,
    ExperimentMetadataResponse,
    TargetMetaData
} from 'src/api/client';
import {ExperimentListItem} from "src/components/types";

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
        async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, errors: string[] }> {
            const response = await getExperiments();
            this.experiments = response;
            const experiments: ExperimentListItem[] = [];
            for (let i = 0; i < response.length; i++) {
                experiments.push({
                    id: response[i].id,
                    name: response[i].name
                });
            }
            return {
                experiments: experiments,
                errors: []
            }
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: [] }> {
            const response = await createExperiment();
            this.experiments.push(response);
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, errors: []
            };
        },
        async deleteExperiment(experimentId: string) {
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
                    experiment_id: experimentId,
                    fasta: targetFile
                };
                const response = await uploadTarget(targetData);
                response.result.map(target => (
                    this.targets.push(target)
                ));
                return response.result;
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
                return await uploadLigand(ligandData);
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
        fetchTargetsForExperiment: async function (experimentId: string) {
            try {
                this.targets = await getTargetsList(experimentId);
            } catch (error) {
                console.error('Error fetching targets:', error);
            }
        },
        async fetchLigandsForTarget(experimentId: string, targetId: string) {
            try {
                return await getLigandsList(experimentId, targetId);
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
        async fetchPocketForTarget(experimentId: string, targetId: string) {
            try {
                const response = await getTargetBindingPocket(experimentId, targetId);
                return response.pocket_ids;
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
