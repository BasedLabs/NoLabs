import {defineStore} from "pinia";
import {
    changeExperimentName,
    deleteExperiment,
    getExperiments,
    getExperiment,
    inference,
    createExperiment
} from "src/features/proteinDesign/api";
import {Experiment, ExperimentListItem, InferenceRequest, ExperimentProperties} from "src/features/proteinDesign/types";
import {Notify} from "quasar";
import {ErrorCodes, ErrorResponse} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";

type ErrorState = {
    error: string
}

const useProteinDesignStore = defineStore("proteinDesign", {
    actions: {
        async inference(request: InferenceRequest): Promise<{
            experiment: Experiment | null,
            error: ErrorState | null
        }> {
            const response = await inference(request.pdbFile, request.contig, request.numberOfDesigns, request.timesteps,
                request.hotspots, request.experimentName, request.experimentId);
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                debugger;
                for(const error of errorResponse.errors){
                    Notify.create({
                        type: "negative",
                        closeBtn: 'Close',
                        message: error
                    });
                }
                return {experiment: null, error: {error: errorResponse.errors[0]}};
            }
            debugger;
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
                error: null
            };
        },
        async getExperiment(experimentId: string): Promise<{
            experiment: Experiment | null,
            error: ErrorState | null
        }> {
            const response = await getExperiment(experimentId);
            const errorResponse = obtainErrorResponse(response);
            if (errorResponse) {
                if (errorResponse.error_code === ErrorCodes.experiment_id_not_found) {
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
                        }, error: null
                    };
                } else {
                    Notify.create({
                        type: "negative",
                        message: errorResponse.errors[0]
                    });
                }

                return {experiment: null, error: {error: errorResponse.errors[0]}};
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
                        numberOfDesigns: response.properties.number_of_desings,
                        timesteps: response.properties.timesteps ?? 50,
                        hotspots: response.properties.hotspots ?? ''
                    }
                }, error: null
            };
        },
        async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, error: ErrorState | null }> {
            const response = await getExperiments();
            const experiments: ExperimentListItem[] = [];
            for (let i = 0; i < response.length; i++) {
                experiments.push({
                    id: response[i].id,
                    name: response[i].name
                });
            }
            return {
                experiments: experiments,
                error: null
            }
        },
        async deleteExperiment(experimentId: string) {
            await deleteExperiment(experimentId);
        },
        async changeExperimentName(experimentId: string, newName: string) {
            await changeExperimentName(experimentId, newName);
        },
        async createExperiment(): Promise<{ experiment: ExperimentListItem | null, error: ErrorState | null }> {
            const response = await createExperiment();
            return {
                experiment: {
                    id: response.id,
                    name: response.name
                }, error: null
            }
        }
    }
});

export default useProteinDesignStore;
