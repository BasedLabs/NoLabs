import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {OpenAPI, LocalisationService, FoldingService} from "src/api/client";
import apiConstants from "src/api/constants";
import {ExperimentListItem} from "src/components/types";
import {Experiment, InferenceRequest} from "src/features/aminoAcid/types";
import {AminoAcid} from "src/features/aminoAcid/folding/types";

OpenAPI.BASE = apiConstants.hostname;

const useFoldingStore = defineStore("folding", {
    actions: {
        async inference(request: InferenceRequest): Promise<{
            experiment: Experiment<AminoAcid> | null,
            errors: string[]
        }> {
            const response = await FoldingService.inferenceApiV1FoldingInferencePost(
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
                            pdbFile: new File([new Blob([aa.pdb_file])], aa.pdb_file_name)
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
            const response = await FoldingService.getExperimentApiV1FoldingGetExperimentGet(experimentId);
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
                    aminoAcids:response.amino_acids.map(aa => {
                        return {
                            name: aa.name,
                            sequence: aa.sequence,
                            pdbFile: new File([new Blob([aa.pdb_file])], aa.pdb_file_name)
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
            const response = await FoldingService.experimentsApiV1FoldingExperimentsGet();
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
            await FoldingService.deleteExperimentApiV1FoldingDeleteExperimentDelete(experimentId);
        },
        async changeExperimentName(experimentId: string, newName: string) {
            await FoldingService.changeExperimentNameApiV1FoldingChangeExperimentNamePost({
                id: experimentId,
                name: newName
            });
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: [] }> {
            const response = await FoldingService.createExperimentApiV1FoldingCreateExperimentGet();
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, errors: []
            }
        }
    }
});

export default useFoldingStore;
