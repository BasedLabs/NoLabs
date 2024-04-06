import {defineStore} from "pinia";
import {Experiment, InferenceRequest} from "src/features/proteinDesign/types";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {OpenAPI, ProteinDesignService} from "src/api/client";
import apiConstants from "src/api/constants";
import {ExperimentListItem} from "src/components/types";

OpenAPI.BASE = apiConstants.hostname;

const useProteinDesignStore = defineStore("proteinDesign", {
    actions: {
        async inference(request: InferenceRequest): Promise<{
            experiment: Experiment | null,
            errors: string[]
        }> {
            const response = await ProteinDesignService.inferenceApiV1ProteinDesignInferencePost(
                {
                    pdb_file: request.pdbFile,
                    contig: request.contig,
                    number_of_designs: request.numberOfDesigns,
                    timesteps: request.timesteps,
                    hotspots: request.hotspots,
                    experiment_name: request.experimentName,
                    experiment_id: request.experimentId
                }
            );
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                for (const error of errorResponse.errors) {
                    Notify.create({
                        type: "negative",
                        closeBtn: 'Close',
                        message: error
                    });
                }
                return {experiment: null, errors: errorResponse.errors};
            }
            return {
                experiment: {
                    id: response.experiment_id,
                    name: response.experiment_name,
                    generatedPdbs: response.pdb_files.map((content, i) => new File([new Blob([content])], `Binder_${i}.pdb`)),
                    properties: {
                        inputPdbFile: request.pdbFile,
                        contig: request.contig,
                        numberOfDesigns: request.numberOfDesigns,
                        timesteps: request.timesteps,
                        hotspots: request.hotspots
                    }
                },
                errors: []
            };
        },
        async getExperiment(experimentId: string): Promise<{
            experiment: Experiment | null,
            errors: string[]
        }> {
            const response = await ProteinDesignService.getExperimentApiV1ProteinDesignExperimentGet(experimentId);
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                if (errorResponse.error_code === ErrorCodes.experiment_not_found) {
                    return {
                        experiment: {
                            id: experimentId,
                            name: "New experiment",
                            generatedPdbs: [],
                            properties: {
                                inputPdbFile: null,
                                contig: '',
                                numberOfDesigns: 2,
                                timesteps: 50,
                                hotspots: ''
                            }
                        }, errors: []
                    };
                } else {
                    Notify.create({
                        type: "negative",
                        message: errorResponse.errors[0]
                    });
                }

                return {experiment: null, errors: errorResponse.errors};
            }

            return {
                experiment: {
                    id: response.experiment_id,
                    name: response.experiment_name,
                    generatedPdbs: response.pdb_files.map((content, i) => new File([new Blob([content])], `Binder_${i}.pdb`)),
                    properties: {
                        inputPdbFile: new File([new Blob([response.properties.pdb_file], {
                            type: 'text/plain'
                        })], response.properties.pdb_file_name),
                        contig: response.properties.contig,
                        numberOfDesigns: response.properties.number_of_designs,
                        timesteps: response.properties.timesteps ?? 50,
                        hotspots: response.properties.hotspots ?? ''
                    }
                }, errors: []
            };
        },
        async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, errors: string[] }> {
            const response = await ProteinDesignService.experimentsApiV1ProteinDesignExperimentsMetadataGet();
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
        async deleteExperiment(experimentId: string) {
            await ProteinDesignService.deleteExperimentApiV1ProteinDesignExperimentDelete(experimentId);
        },
        async changeExperimentName(experimentId: string, newName: string) {
            await ProteinDesignService.changeExperimentNameApiV1ProteinDesignChangeExperimentNamePost({
                id: experimentId,
                name: newName
            });
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: [] }> {
            const response = await ProteinDesignService.createExperimentApiV1ProteinDesignCreateExperimentGet();
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, errors: []
            }
        }
    }
});

export default useProteinDesignStore;
