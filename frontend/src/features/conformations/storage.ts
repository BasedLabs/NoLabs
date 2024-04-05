import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes, ErrorResponse} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {Experiment, InferenceRequest} from "src/features/conformations/types";
import {ConformationsService, IntegratorsRequest, OpenAPI} from "src/api/client";
import apiConstants from "src/api/constants";
import {ExperimentListItem} from "src/components/types";


OpenAPI.BASE = apiConstants.hostname;

const useConformationsStore = defineStore("conformations", {
    actions: {
        async inference(request: InferenceRequest): Promise<{
            experiment: Experiment | null,
            errors: string[]
        }> {
            const response = await ConformationsService.inferenceApiV1ConformationsInferencePost(
                {
                    pdb_file: request.pdbFile,
                    experiment_name: request.experimentName,
                    experiment_id: request.experimentId,
                    total_frames: request.totalFrames,
                    temperature_k: request.temperatureK,
                    take_frame_every: request.takeFrameEvery,
                    step_size: request.stepSize,
                    replace_non_standard_residues: request.replaceNonStandardResidues,
                    add_missing_atoms: request.addMissingAtoms,
                    add_missing_hydrogens: request.addMissingHydrogens,
                    friction_coeff: request.frictionCoeff,
                    ignore_missing_atoms: request.ignoreMissingAtoms,
                    integrator: request.integrator,
                }
            )
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
                    pdbContent: new File([new Blob([response.pdb_content!])], request.pdbFile.name),
                    timeline: response.timeline.map(x => {
                        return {
                            message: x.message,
                            error: x.error,
                            createdAt: x.created_at
                        };
                    }),
                    properties: {
                        pdbFile: request.pdbFile,
                        totalFrames: request.totalFrames,
                        temperatureK: request.temperatureK,
                        takeFrameEvery: request.takeFrameEvery,
                        stepSize: request.stepSize,
                        replaceNonStandardResidues: request.replaceNonStandardResidues,
                        addMissingAtoms: request.addMissingAtoms,
                        addMissingHydrogens: request.addMissingHydrogens,
                        frictionCoeff: request.frictionCoeff,
                        ignoreMissingAtoms: request.ignoreMissingAtoms,
                        integrator: request.integrator,
                    }
                },
                errors: []
            };
        },
        async getExperiment(experimentId: string): Promise<{
            experiment: Experiment | null,
            errors: string[]
        }> {
            const response = await ConformationsService.getExperimentApiV1ConformationsGetExperimentGet(experimentId);
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                if (errorResponse.error_code === ErrorCodes.experiment_not_found) {
                    return {
                        experiment: {
                            id: experimentId,
                            name: 'New experiment',
                            pdbContent: null,
                            timeline: [],
                            properties: {
                                pdbFile: null,
                                totalFrames: 10000,
                                temperatureK: 273.15,
                                takeFrameEvery: 1000,
                                stepSize: 0.002,
                                replaceNonStandardResidues: false,
                                addMissingAtoms: false,
                                addMissingHydrogens: true,
                                frictionCoeff: 1.0,
                                ignoreMissingAtoms: false,
                                integrator: IntegratorsRequest.LANGEVIN_INTEGATOR
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
                    pdbContent: new File([new Blob([response.pdb_file!])], response.properties.pdb_file_name),
                    timeline: response.timeline.map(x => {
                        return {
                            message: x.message,
                            error: x.error,
                            createdAt: x.created_at
                        };
                    }),
                    properties: {
                        pdbFile: new File([new Blob([response.properties.pdb_file!])], response.properties.pdb_file_name),
                        totalFrames: response.properties.total_frames,
                        temperatureK: response.properties.temperature_k,
                        takeFrameEvery: response.properties.take_frame_every,
                        stepSize: response.properties.step_size,
                        replaceNonStandardResidues: response.properties.replace_non_standard_residues,
                        addMissingAtoms: response.properties.add_missing_atoms,
                        addMissingHydrogens: response.properties.add_missing_hydrogens,
                        frictionCoeff: response.properties.friction_coeff,
                        ignoreMissingAtoms: response.properties.ignore_missing_atoms,
                        integrator: response.properties.integrator!
                    }
                }, errors: []
            };
        },
        async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, errors: string[] }> {
            const response = await ConformationsService.experimentsApiV1ConformationsExperimentsGet();
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
            await ConformationsService.deleteExperimentApiV1ConformationsDeleteExperimentDelete(experimentId);
        },
        async changeExperimentName(experimentId: string, newName: string) {
            await ConformationsService.changeExperimentNameApiV1ConformationsChangeExperimentNamePost({
                id: experimentId,
                name: newName
            });
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: string[] }> {
            const response = await ConformationsService.createExperimentApiV1ConformationsCreateExperimentGet();
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, errors: []
            }
        }
    }
});

export default useConformationsStore;
