import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {OpenAPI, SolubilityService} from "src/api/client";
import apiConstants from "src/api/constants";
import {ExperimentListItem} from "src/components/types";
import {Experiment, InferenceRequest} from "src/features/aminoAcid/types";
import {AminoAcid} from "src/features/aminoAcid/solubility/types";

OpenAPI.BASE = apiConstants.hostname;

const useSolubilityStore = defineStore("solubility", {
    actions: {
        async inference(request: InferenceRequest): Promise<{
            experiment: Experiment<AminoAcid> | null,
            errors: string[]
        }> {
            const response = await SolubilityService.inferenceApiV1SolubilityInferencePost(
                {
                    experiment_id: request.experimentId,
                    experiment_name: request.experimentName,
                    amino_acid_sequence: request.aminoAcidSequence,
                    fastas: request.fastas
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
                    aminoAcids: response.amino_acids.map(aa => {
                        return {
                            name: aa.name,
                            sequence: aa.sequence,
                            soluble_probability: aa.soluble_probability
                        };
                    }),
                    properties: {
                        aminoAcidSequence: request.aminoAcidSequence,
                        fastas: request.fastas
                    }
                },
                errors: []
            };
        },
        async getExperiment(experimentId: string): Promise<{
            experiment: Experiment<AminoAcid> | null,
            errors: string[]
        }> {
            const response = await SolubilityService.getExperimentApiV1SolubilityGetExperimentGet(experimentId);
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                if (errorResponse.error_code === ErrorCodes.experiment_not_found) {
                    return {
                        experiment: {
                            id: experimentId,
                            name: "New experiment",
                            aminoAcids: [],
                            properties: {
                                aminoAcidSequence: null,
                                fastas: []
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
                    aminoAcids: response.amino_acids.map(aa => {
                        return {
                            name: aa.name,
                            sequence: aa.sequence,
                            soluble_probability: aa.soluble_probability
                        };
                    }),
                    properties: {
                        aminoAcidSequence: response.properties.amino_acid_sequence as (string | undefined),
                        fastas: response.properties.fastas.map(f =>
                            new File([new Blob([f.content])], f.filename)
                        )
                    }
                }, errors: []
            };
        },
        async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, errors: string[] }> {
            const response = await SolubilityService.experimentsApiV1SolubilityExperimentsGet();
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
            await SolubilityService.deleteExperimentApiV1SolubilityDeleteExperimentDelete(experimentId);
        },
        async changeExperimentName(experimentId: string, newName: string) {
            await SolubilityService.changeExperimentNameApiV1SolubilityChangeExperimentNamePost({
                id: experimentId,
                name: newName
            });
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: [] }> {
            const response = await SolubilityService.createExperimentApiV1SolubilityCreateExperimentGet();
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, errors: []
            }
        }
    }
});

export default useSolubilityStore;
