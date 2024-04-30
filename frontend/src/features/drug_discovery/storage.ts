import {defineStore} from 'pinia';
import {
  changeExperimentNameApi,
  changeTargetNameApi,
  checkDiffDockJobIsRunningApi,
  checkDiffDockServiceHealthApi,
  checkDockingResultAvailableApi,
  checkEsmFoldJobIsRunningApi,
  checkEsmFoldLightJobIsRunningApi,
  checkEsmFoldLightServiceHealthApi,
  checkEsmFoldServiceHealthApi,
  checkFoldingDataAvailableApi,
  checkMsaDataAvailableApi,
  checkMsaJobIsRunningApi,
  checkMsaServiceHealthApi,
  checkP2RankJobIsRunningApi,
  checkP2RankServiceHealthApi,
  checkPocketDataAvailableApi, checkRosettaFoldJobIsRunningApi, checkRosettaFoldServiceHealthApi,
  checkUmolJobIsRunningApi,
  checkUmolServiceHealthApi,
  createExperimentApi,
  deleteDockingJobApi,
  deleteExperimentApi,
  deleteLigandFromExperimentApi,
  deleteLigandFromTargetApi,
  deleteTargetApi,
  getAllDockingJobsListApi,
  getAllDockingResultsListApi,
  getDiffDockDockingJobResultDataApi,
  getDiffDockLigandSdfApi,
  getDiffDockParamsApi,
  getDockingJobsListForTargetLigandApi,
  getDockingResultsListForTargetLigandApi,
  getExperimentMetadataApi,
  getExperimentsApi,
  getJobPocketDataApi,
  getLigandDataForExperimentApi,
  getLigandDataForTargetApi,
  getLigandMetaDataForExperimentApi,
  getLigandMetaDataForTargetApi,
  getLigandsListForExperimentApi,
  getLigandsListForTargetApi,
  getTargetBindingPocketApi,
  getTargetDataApi,
  getTargetMetaDataApi,
  getTargetsListApi,
  getUmolDockingJobResultDataApi,
  predictBindingPocketApi,
  predictFoldingApi,
  predictLightFoldingApi, predictRosettaFoldApi,
  registerDockingJobApi,
  runDiffDockDockingJobApi,
  runUmolDockingJobApi,
  setTargetBindingPocketApi,
  updateDiffDockParamsApi,
  updateDockingParamsApi,
  uploadLigandToExperimentApi,
  uploadLigandToTargetApi,
  uploadTargetApi
} from 'src/features/drug_discovery/api';

import {
  Body_upload_ligand_to_experiment_api_v1_drug_discovery_upload_ligand_to_experiment_post,
  Body_upload_ligand_to_target_api_v1_drug_discovery_upload_ligand_to_target_post,
  Body_upload_target_api_v1_drug_discovery_upload_target_post,
  ExperimentMetadataResponse,
} from 'src/api/client';
import {ExperimentListItem} from "src/components/types";
import {obtainErrorResponse} from "../../api/errorWrapper";
import {Notify} from "quasar";
import {ExtendedLigandMetaData, ExtendedTargetMetaData} from "./components/targets/types";

export const useDrugDiscoveryStore = defineStore('drugDiscovery', {

    state: () => ({
        experiments: [] as ExperimentMetadataResponse[],
        currentExperiment: null as ExperimentMetadataResponse | null,
        targets: [] as ExtendedTargetMetaData[],
        loneLigands: [] as ExtendedLigandMetaData[],
        currentLigand: null,
        foldingMethods: {
            esmfoldLight: 'light',
            esmfold: 'esmfold',
            rosettafold: 'rosettafold'
        }
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
        async uploadTargetToExperiment(experimentId: string, targetFile: File, metaData?: Record<string, string>) {
            try {
                const targetData: Body_upload_target_api_v1_drug_discovery_upload_target_post = {
                    experiment_id: experimentId,
                    fasta: targetFile,
                    ...(metaData && { metadata: JSON.stringify(metaData) })
                };
                const response = await uploadTargetApi(targetData);

                let uploadedTargets = [] as ExtendedTargetMetaData[];
                response.result.map(target => (
                  uploadedTargets.push({...target,
                      data: undefined,
                      loadingLigands: false,
                      ligands: [],
                      loadingTargetData: false})
                ));
                return uploadedTargets;
            } catch (error) {
                console.error('Error uploading target:', error);
            }
        },
        async deleteTargetFromExperiment(experimentId: string, targetId: string) {
            try {
                await deleteTargetApi(experimentId, targetId);
                this.targets = this.targets.filter(target => target.target_id !== targetId);
            } catch (error) {
                console.error('Error deleting target:', error);
            }
        },
        async uploadLigandToTarget(experimentId: string, targetId: string, sdfFile: File) {
            try {
                const ligandData: Body_upload_ligand_to_target_api_v1_drug_discovery_upload_ligand_to_target_post = {
                    experiment_id: experimentId,
                    target_id: targetId,
                    sdf_file: sdfFile
                };
                const response = await uploadLigandToTargetApi(ligandData);
                return response.ligand_meta_data as ExtendedLigandMetaData;
            } catch (error) {
                console.error('Error uploading ligand:', error);
            }
        },
        async uploadLigandToExperiment(experimentId: string, sdfFile: File, metaData?: Record<string, string>) {
          try {
            const ligandData: Body_upload_ligand_to_experiment_api_v1_drug_discovery_upload_ligand_to_experiment_post = {
              experiment_id: experimentId,
              sdf_file: sdfFile,
              ...(metaData && { metadata: JSON.stringify(metaData) })
            };
            const response = await uploadLigandToExperimentApi(ligandData);
            let newLigand = {
              ...response.ligand_meta_data,
              jobs: [],
              loadingLigandData: false
            } as ExtendedLigandMetaData;
            newLigand.data = await this.fetchLigandDataForExperiment(experimentId, response.ligand_meta_data.ligand_id);
            return newLigand;
          } catch (error) {
            console.error('Error uploading ligand:', error);
          }
        },
        async deleteLigandFromTarget(experimentId: string, targetId: string, ligandId: string) {
            try {
                await deleteLigandFromTargetApi(experimentId, targetId, ligandId);
                // Update the ligands list for the specific target
            } catch (error) {
                console.error('Error deleting ligand:', error);
            }
        },
        async deleteLigandFromExperiment(experimentId: string, ligandId: string) {
          try {
            await deleteLigandFromExperimentApi(experimentId, ligandId);
            this.loneLigands = this.loneLigands.filter(ligand => ligand.ligand_id !== ligandId);
          } catch (error) {
            console.error('Error deleting ligand:', error);
          }
        },
        fetchTargetsForExperiment: async function (experimentId: string) {
            try {
                return (await getTargetsListApi(experimentId)).map(x => {return { ...x,
                  data: undefined,
                  loadingLigands: false,
                  ligands: [],
                  loadingTargetData: false}}) as ExtendedTargetMetaData[];
            } catch (error) {
                console.error('Error fetching targets:', error);
            }
        },
        async fetchLigandsForTarget(experimentId: string, targetId: string) {
            try {
                return await getLigandsListForTargetApi(experimentId, targetId);
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchLigandsForExperiment(experimentId: string) {
          try {
            return (await getLigandsListForExperimentApi(experimentId)).map(x => {return { ...x,
              jobs: [],
              loadingLigandData: false}}) as ExtendedLigandMetaData[];
          } catch (error) {
            console.error('Error fetching ligands:', error);
          }
        },
        async fetchTargetMetaData(experimentId: string, targetId: string) {
            try {
                return await getTargetMetaDataApi(experimentId, targetId);
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async changeTargetName(experimentId: string, targetId: string, targetName: string) {
            try {
                return await changeTargetNameApi(experimentId, targetId, targetName);
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchTargetData(experimentId: string, targetId: string) {
            try {
                const response = await getTargetDataApi(experimentId, targetId);

                return {proteinSequence: response.protein_sequence,
                    pdbContents: response.protein_pdb };
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async predictEsmLightFoldingForTarget(experimentId: string, targetId: string) {
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
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async predictEsmFoldingForTarget(experimentId: string, targetId: string) {
            try {
                const response = await predictFoldingApi(experimentId, targetId);
                return response.pdb_content;
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchPocketForTarget(experimentId: string, targetId: string) {
            try {
                const response = await getTargetBindingPocketApi(experimentId, targetId);
                return response.pocket_ids;
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
       async predictRoseTTAFold(experimentId: string, targetId: string) {
            try {
                const response = await predictRosettaFoldApi(experimentId, targetId);
                return response.pdb_content;
                // Optionally set the currentTarget to the selected one
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async setPocketForTarget(experimentId: string, targetId: string, pocketIds: Array<number>) {
            try {
                await setTargetBindingPocketApi(experimentId, targetId, pocketIds);
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async predictPocketForTarget(experimentId: string, targetId: string) {
            try {
                const response = await predictBindingPocketApi(experimentId, targetId);
                return response.pocket_ids;
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchLigandMetaDataForTarget(experimentId: string, targetId: string, ligandId: string) {
            try {
                return await getLigandMetaDataForTargetApi(experimentId, targetId, ligandId);
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchLigandMetaDataForExperiment(experimentId: string, ligandId: string) {
          try {
            return await getLigandMetaDataForExperimentApi(experimentId, ligandId);
          } catch (error) {
            console.error('Error fetching ligands:', error);
          }
        },
        async fetchLigandDataForTarget(experimentId: string, targetId: string, ligandId: string) {
            try {
                const response = await getLigandDataForTargetApi(experimentId, targetId, ligandId);
                return {
                    ligandId: ligandId,
                    ligandName: response.ligand_name,
                    ligandSdf: response.ligand_sdf,
                    ligandSmiles: response.ligand_smiles
                };
            } catch (error) {
                console.error('Error fetching ligands:', error);
            }
        },
        async fetchLigandDataForExperiment(experimentId: string, ligandId: string) {
          try {
            const response = await getLigandDataForExperimentApi(experimentId, ligandId);
            return {
                sdf_file: response.ligand_sdf,
                smiles: response.ligand_smiles
              };
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
        async updateDockingParams(experimentId: string, targetId: string, ligandId: string, jobId: string, foldingMethod: string, dockingMethod: string) {
            try {
                return await updateDockingParamsApi(experimentId, targetId, ligandId, jobId, foldingMethod, dockingMethod);
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

        async checkRosettafoldJobIsRunning(jobId: string) {
            try {
                return await checkRosettaFoldJobIsRunningApi(jobId);
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
        async checkRosettaFoldServiceHealth() {
            try {
                return await checkRosettaFoldServiceHealthApi();
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
