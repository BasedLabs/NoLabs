import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {GeneOntologyService, OpenAPI, RunGeneOntologyResponse, SolubilityService} from "src/api/client";
import apiConstants from "src/api/constants";
import {ExperimentListItem} from "src/components/types";
import {Experiment, InferenceRequest} from "src/features/aminoAcid/types";
import {AminoAcid} from "src/features/aminoAcid/geneOntology/types";

OpenAPI.BASE = apiConstants.hostname;

const useGeneOntologyStore = defineStore("geneOntology", {
    actions: {
        _mapResponse(response: RunGeneOntologyResponse, aminoAcidSequence: string | undefined | null, fastas: Array<File>): Experiment<AminoAcid> {
            return {
                id: response.experiment_id,
                name: response.experiment_name,
                aminoAcids: response.amino_acids.map(aa => {
                    return {
                        name: aa.name,
                        sequence: aa.sequence,
                        go: aa.go
                    }
                }),
                properties: {
                    aminoAcidSequence: aminoAcidSequence,
                    fastas: fastas
                }
            }
        },
        async inference(request: InferenceRequest): Promise<{
            experiment: Experiment<AminoAcid> | null,
            errors: string[]
        }> {
            const response = await GeneOntologyService.inferenceApiV1GeneOntologyInferencePost(
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
                experiment: this._mapResponse(response, request.aminoAcidSequence, request.fastas),
                errors: []
            };
        },
        async getExperiment(experimentId: string): Promise<{
            experiment: Experiment<AminoAcid> | null,
            errors: string[]
        }> {
            const response = await GeneOntologyService.getExperimentApiV1GeneOntologyGetExperimentGet(experimentId);
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                if (errorResponse.error_code === ErrorCodes.experiment_not_found) {
                    return {
                        experiment: {
                            id: experimentId,
                            name: 'New experiment',
                            aminoAcids: [],
                            properties: {
                                aminoAcidSequence: undefined,
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
                experiment: this._mapResponse(response,
                    response.properties.amino_acid_sequence,
                    response.properties.fastas.map(f =>
                        new File([new Blob([f.content])], f.filename)
                    )),
                errors: []
            };
        },
        async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, errors: string[] }> {
            const response = await GeneOntologyService.experimentsApiV1GeneOntologyExperimentsGet();
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
            await GeneOntologyService.deleteExperimentApiV1GeneOntologyDeleteExperimentDelete(experimentId);
        },
        async changeExperimentName(experimentId: string, newName: string) {
            await GeneOntologyService.changeExperimentNameApiV1GeneOntologyChangeExperimentNamePost({
                id: experimentId,
                name: newName
            });
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: [] }> {
            const response = await GeneOntologyService.createExperimentApiV1GeneOntologyCreateExperimentGet();
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, errors: []
            }
        }
    }
});

export default useGeneOntologyStore;
