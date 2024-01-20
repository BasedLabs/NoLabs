import {
    RunProteinDesignResponse,
    ProteinDesignService,
    CancelablePromise,
    nolabs__api_models__protein_design__GetExperimentResponse,
    nolabs__api_models__protein_design__ExperimentMetadataResponse
} from '../../api/client';

function inference(
    pdbFile: Blob,
    contig: string,
    numberOfDesigns: number,
    timesteps: number,
    hotspots: string,
    experimentName: string,
    experimentId?: string
): CancelablePromise<RunProteinDesignResponse> {
    const formData = {
        experiment_name: experimentName,
        experiment_id: experimentId,
        pdb_file: pdbFile,
        contig: contig,
        number_of_designs: numberOfDesigns,
        timesteps: timesteps,
        hotspots: hotspots
    }
    return ProteinDesignService.inferenceApiV1ProteinDesignInferencePost(formData)
}

function getExperiment(experimentId: string): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
    return ProteinDesignService.getExperimentApiV1ProteinDesignGetExperimentGet(experimentId);
}

function getExperiments(): CancelablePromise<Array<nolabs__api_models__protein_design__ExperimentMetadataResponse>> {
    return ProteinDesignService.experimentsApiV1ProteinDesignExperimentsGet();
}

function deleteExperiment(experimentId: string): CancelablePromise<any> {
    return ProteinDesignService.deleteExperimentApiV1ProteinDesignDeleteExperimentDelete(experimentId);
}

function changeExperimentName(experimentId: string, experimentName: string): CancelablePromise<any> {
    return ProteinDesignService.changeExperimentNameApiV1ProteinDesignChangeExperimentNamePost({
        id: experimentId,
        name: experimentName
    })
}

export { inference, getExperiment, getExperiments, deleteExperiment, changeExperimentName };