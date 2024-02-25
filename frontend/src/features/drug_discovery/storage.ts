import {defineStore} from 'pinia';
import {
  changeExperimentNameApi,
  changeTargetNameApi,
  checkDockingResultAvailableApi,
  checkFoldingDataAvailableApi,
  checkEsmFoldJobIsRunningApi,
  checkEsmFoldLightJobIsRunningApi,
  checkEsmFoldServiceHealthApi,
  checkEsmFoldLightServiceHealthApi,
  checkMsaDataAvailableApi,
  checkMsaJobIsRunningApi,
  checkMsaServiceHealthApi,
  checkP2RankJobIsRunningApi,
  checkP2RankServiceHealthApi,
  checkPocketDataAvailableApi,
  checkUmolJobIsRunningApi,
  checkUmolServiceHealthApi,
  createExperimentApi,
  deleteDockingJobApi,
  deleteExperimentApi,
  deleteLigandApi,
  deleteTargetApi,
  getAllDockingJobsListApi,
  getAllDockingResultsListApi,
  getUmolDockingJobResultDataApi,
  getDockingJobsListForTargetLigandApi,
  getDockingResultsListForTargetLigandApi,
  getExperimentMetadataApi,
  getExperimentsApi,
  getJobPocketDataApi,
  getLigandDataApi,
  getLigandMetaDataApi,
  getLigandsListApi,
  getTargetBindingPocketApi,
  getTargetDataApi,
  getTargetMetaDataApi,
  getTargetsListApi,
  predictBindingPocketApi,
  predictFoldingApi,
  predictLightFoldingApi,
  registerDockingJobApi,
  runUmolDockingJobApi,
  setTargetBindingPocketApi,
  uploadLigandApi,
  uploadTargetApi,
  checkDiffDockServiceHealthApi,
  checkDiffDockJobIsRunningApi,
  runDiffDockDockingJobApi,
  getDiffDockDockingJobResultDataApi,
  getDiffDockLigandSdfApi,
  getDiffDockParamsApi,
  updateDiffDockParamsApi,
  updateDockingParamsApi
} from 'src/features/drug_discovery/api';

import {
    Body_upload_ligand_api_v1_drug_discovery_upload_ligand_post,
    Body_upload_target_api_v1_drug_discovery_upload_target_post,
    ExperimentMetadataResponse,
    TargetMetaData
} from 'src/api/client';
import {ExperimentListItem} from "src/components/types";
import {obtainErrorResponse} from "../../api/errorWrapper";
import {Notify} from "quasar";

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
            const response = await getExperimentsApi();
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
            const response = await createExperimentApi();
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
                await deleteExperimentApi(experimentId);
                this.experiments = this.experiments.filter(exp => exp.id !== experimentId);
            } catch (error) {
                console.error('Error deleting experiment:', error);
            }
        },
        async changeExperimentName(experimentId: string, experimentName: string) {
            try {
                return await changeExperimentNameApi(experimentId, experimentName);
            } catch (error) {
                console.error('Error deleting experiment:', error);
            }
        },
        async getExperimentMetaData(experimentId: string) {
            try {
                return await getExperimentMetadataApi(experimentId);
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
                const response = await uploadTargetApi(targetData);
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
                await deleteTargetApi(experimentId, targetId);
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
                const response = await uploadLigandApi(ligandData);
                return {
                    ligand_id: response.ligand_meta_data.ligand_id,
                    ligand_name: response.ligand_meta_data.ligand_name,
                    jobs: []
                }
            } catch (error) {
                console.error('Error uploading ligand:', error);
            }
        },
        async deleteLigandFromTarget(experimentId: string, targetId: string, ligandId: string) {
            try {
                await deleteLigandApi(experimentId, targetId, ligandId);
                // Update the ligands list for the specific target
            } catch (error) {
                console.error('Error deleting ligand:', error);
            }
        },
        fetchTargetsForExperiment: async function (experimentId: string) {
            try {
                this.targets = (await getTargetsListApi(experimentId)).map(x => {return {...x, ligands: []}});
            } catch (error) {
                console.error('Error fetching targets:', error);
            }
        },
        async fetchLigandsForTarget(experimentId: string, targetId: string) {
            try {
                return await getLigandsListApi(experimentId, targetId);
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchTargetMetaData(experimentId: string, targetId: string) {
            try {
                return await getTargetMetaDataApi(experimentId, targetId);
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async changeTargetName(experimentId: string, targetId: string, targetName: string) {
            try {
                return await changeTargetNameApi(experimentId, targetId, targetName);
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchTargetData(experimentId: string, targetId: string) {
            try {
                const response = await getTargetDataApi(experimentId, targetId);
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
        async predictLightFoldingForTarget(experimentId: string, targetId: string) {
            try {
                const response = await predictLightFoldingApi(experimentId, targetId);
                const errorResponse = obtainErrorResponse(response);
                if (errorResponse) {
                    for (const error of errorResponse.errors) {
                        Notify.create({
                            type: "negative",
                            closeBtn: 'Close',
                            message: error
                        });
                    }
                }
                return response.pdb_content;
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async predictFoldingForTarget(experimentId: string, targetId: string) {
            try {
                const response = await predictFoldingApi(experimentId, targetId);
                return response.pdb_content;
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchPocketForTarget(experimentId: string, targetId: string, folding_method: string) {
            try {
                const response = await getTargetBindingPocketApi(experimentId, targetId, folding_method);
                return response.pocket_ids;
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async setPocketForTarget(experimentId: string, targetId: string, pocketIds: Array<number>) {
            try {
                await setTargetBindingPocketApi(experimentId, targetId, pocketIds);
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async predictPocketForTarget(experimentId: string, targetId: string, folding_method: string) {
            try {
                const response = await predictBindingPocketApi(experimentId, targetId, folding_method);
                return response.pocket_ids;
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchLigandMetaData(experimentId: string, targetId: string, ligandId: string) {
            try {
                return await getLigandMetaDataApi(experimentId, targetId, ligandId);
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchLigandData(experimentId: string, targetId: string, ligandId: string) {
            try {
                const response = await getLigandDataApi(experimentId, targetId, ligandId);
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
        async registerDockingJob(experimentId: string, targetId: string, ligandId: string, foldingMethod: string) {
            try {
                return await registerDockingJobApi(experimentId, targetId, ligandId, foldingMethod);
            } catch (error) {
                console.error("Error registering docking job:", error);
                return null;
            }
        },

        async runUmolDockingJob(experimentId: string, targetId: string, ligandId: string, jobId: string) {
            try {
                return await runUmolDockingJobApi(experimentId, targetId, ligandId, jobId);
            } catch (error) {
                console.error("Error running docking job:", error);
                return null;
            }
        },

        async runDiffDockDockingJob(experimentId: string, targetId: string, ligandId: string, jobId: string) {
          try {
            return await runDiffDockDockingJobApi(experimentId, targetId, ligandId, jobId);
          } catch (error) {
            console.error("Error running docking job:", error);
            return null;
          }
        },

        async getUmolDockingJobResultData(experimentId: string, targetId: string, ligandId: string, jobId: string) {
            try {
                return await getUmolDockingJobResultDataApi(experimentId, targetId, ligandId, jobId);
            } catch (error) {
                console.error("Error getting docking job result data:", error);
                return null;
            }
        },

        async getDiffDockDockingJobResultData(experimentId: string, targetId: string, ligandId: string, jobId: string) {
          try {
            return await getDiffDockDockingJobResultDataApi(experimentId, targetId, ligandId, jobId);
          } catch (error) {
            console.error("Error getting docking job result data:", error);
            return null;
          }
        },

        async getDiffDockParams(experimentId: string, targetId: string, ligandId: string, jobId: string) {
          try {
            return await getDiffDockParamsApi(experimentId, targetId, ligandId, jobId);
          } catch (error) {
            console.error("Error getting docking job result data:", error);
            return null;
          }
        },

        async updateDockingParams(experimentId: string, targetId: string, ligandId: string, jobId: string, folding_method: string, docking_method: string) {
          try {
            return await updateDockingParamsApi(experimentId, targetId, ligandId, jobId, folding_method, docking_method);
          } catch (error) {
            console.error("Error getting docking job result data:", error);
            return null;
          }
        },

        async updateDiffDockParams(experimentId: string, targetId: string, ligandId: string, jobId: string, samples_per_complex: number) {
          try {
            return await updateDiffDockParamsApi(experimentId, targetId, ligandId, jobId, samples_per_complex);
          } catch (error) {
            console.error("Error getting docking job result data:", error);
            return null;
          }
        },

        async getDiffDockLigandSdf(experimentId: string, targetId: string, ligandId: string, jobId: string, ligandFileName: string) {
          try {
            return await getDiffDockLigandSdfApi(experimentId, targetId, ligandId, jobId, ligandFileName);
          } catch (error) {
            console.error("Error getting docking job result data:", error);
            return null;
          }
        },

        async checkDockingResultAvailable(experimentId: string, targetId: string, ligandId: string, jobId: string) {
            try {
                return await checkDockingResultAvailableApi(experimentId, targetId, ligandId, jobId);
            } catch (error) {
                console.error("Error checking docking result availability:", error);
                return null;
            }
        },

        async checkMsaDataAvailable(experimentId: string, targetId: string) {
            try {
                return await checkMsaDataAvailableApi(experimentId, targetId);
            } catch (error) {
                console.error("Error checking MSA data availability:", error);
                return null;
            }
        },

        async checkMsaJobIsRunning(jobId: string) {
            try {
                return await checkMsaJobIsRunningApi(jobId);
            } catch (error) {
                console.error("Error checking if MSA job is running:", error);
                return null;
            }
        },

        async checkFoldingDataAvailable(experimentId: string, targetId: string, foldingMethod: string) {
            try {
                return await checkFoldingDataAvailableApi(experimentId, targetId, foldingMethod);
            } catch (error) {
                console.error("Error checking folding data availability:", error);
                return null;
            }
        },

        async checkEsmFoldJobIsRunning(jobId: string) {
            try {
                return await checkEsmFoldJobIsRunningApi(jobId);
            } catch (error) {
                console.error("Error checking if folding job is running:", error);
                return null;
            }
        },

        async checkEsmFoldLightJobIsRunning(jobId: string) {
          try {
            return await checkEsmFoldLightJobIsRunningApi(jobId);
          } catch (error) {
            console.error("Error checking if folding job is running:", error);
            return null;
          }
        },

        async checkPocketDataAvailable(experimentId: string, targetId: string, ligandId: string, jobId: string) {
            try {
                return await checkPocketDataAvailableApi(experimentId, targetId, ligandId, jobId);
            } catch (error) {
                console.error("Error checking pocket data availability:", error);
                return null;
            }
        },

        async checkP2RankJobIsRunning(jobId: string) {
            try {
                return await checkP2RankJobIsRunningApi(jobId);
            } catch (error) {
                console.error("Error checking if P2Rank job is running:", error);
                return null;
            }
        },
        async checkUmolJobIsRunning(jobId: string) {
            try {
                return await checkUmolJobIsRunningApi(jobId);
            } catch (error) {
                console.error("Error checking if Umol job is running:", error);
                return null;
            }
        },
        async checkDiffDockJobIsRunning(jobId: string) {
          try {
            return await checkDiffDockJobIsRunningApi(jobId);
          } catch (error) {
            console.error("Error checking if Umol job is running:", error);
            return null;
          }
        },
        async getAllDockingResultsList(experimentId: string) {
            try {
                return await getAllDockingResultsListApi(experimentId);
            } catch (error) {
                console.error("Error getting all results list:", error);
                return null;
            }
        },
        async getAllDockingJobsList(experimentId: string) {
          try {
            return await getAllDockingJobsListApi(experimentId);
          } catch (error) {
            console.error("Error getting all results list:", error);
            return null;
          }
        },
        async getDockingResultsListForTargetLigand(experimentId: string, targetId: string, ligandId: string) {
            try {
                return await getDockingResultsListForTargetLigandApi(experimentId, targetId, ligandId);
            } catch (error) {
                console.error("Error getting results list for ligand-target:", error);
                return null;
            }
        },
        async getDockingJobsListForTargetLigand(experimentId: string, targetId: string, ligandId: string) {
          try {
            return await getDockingJobsListForTargetLigandApi(experimentId, targetId, ligandId);
          } catch (error) {
            console.error("Error getting results list for ligand-target:", error);
            return null;
          }
        },
        async deleteDockingJob(experimentId: string, targetId: string, ligandId: string, jobId: string) {
            try {
                return await deleteDockingJobApi(experimentId, targetId, ligandId, jobId);
            } catch (error) {
                console.error("Error deleting docking job:", error);
                return null;
            }
        },
        async getJobPocketIds(experimentId: string, targetId: string, ligandId: string, jobId: string) {
            try {
                return await getJobPocketDataApi(experimentId, targetId, ligandId, jobId);
            } catch (error) {
                console.error("Error getting result pocket Ids:", error);
                return null;
            }
        },
        async checkMsaServiceHealth() {
            try {
                return await checkMsaServiceHealthApi();
            } catch (error) {
                console.error("Error checking if msa service is healthy:", error);
                return null;
            }
        },
        async checkP2RankServiceHealth() {
            try {
                return await checkP2RankServiceHealthApi();
            } catch (error) {
                console.error("Error checking if p2rank service is healthy:", error);
                return null;
            }
        },
        async checkEsmFoldServiceHealth() {
            try {
                return await checkEsmFoldServiceHealthApi();
            } catch (error) {
                console.error("Error checking if folding service is healthy:", error);
                return null;
            }
        },
        async checkEsmFoldLightServiceHealth() {
          try {
            return await checkEsmFoldLightServiceHealthApi();
          } catch (error) {
            console.error("Error checking if folding service is healthy:", error);
            return null;
          }
        },
        async checkUmolServiceHealth() {
            try {
                return await checkUmolServiceHealthApi();
            } catch (error) {
                console.error("Error checking if Umol service is healthy:", error);
                return null;
            }
        },
        async checkDiffDockServiceHealth() {
          try {
            return await checkDiffDockServiceHealthApi();
          } catch (error) {
            console.error("Error checking if Umol service is healthy:", error);
            return null;
          }
        },
    }
});
