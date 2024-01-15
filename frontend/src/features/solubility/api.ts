import {
    RunSolubilityResponse,
    SolubilityService,
    CancelablePromise,
    nolabs__api_models__solubility__GetExperimentResponse,
    ExperimentMetadataResponse
} from 'api/client';

async function inference(
    fastas: Blob[],
    experimentName?: string,
    experimentId?: string,
    aminoAcidSequence?: string): CancelablePromise<RunSolubilityResponse> {
    const formData = {
        experiment_name: experimentName,
        experiment_id: experimentId,
        amino_acid_sequence: aminoAcidSequence,
        fastas: fastas,
    }
    return await SolubilityService.inferenceApiV1SolubilityInferencePost(formData)
}

async function getExperiment(experimentId: string): CancelablePromise<nolabs__api_models__solubility__GetExperimentResponse> {
    return await SolubilityService.getExperimentApiV1SolubilityGetExperimentGet(experimentId);
}

async function experiments(): CancelablePromise<Array<ExperimentMetadataResponse>> {
    return await SolubilityService.experimentsApiV1SolubilityExperimentsGet();
}

async function deleteExperiment(experimentId: string): CancelablePromise<any> {
    return await SolubilityService.deleteExperimentApiV1SolubilityDeleteExperimentDelete(experimentId);
}

async function changeExperimentName(experimentId: string, experimentName: string): CancelablePromise<any> {
    return await SolubilityService.changeExperimentNameApiV1SolubilityChangeExperimentNamePost({
        id: experimentId,
        name: experimentName
    })
}

export { inference, getExperiment, experiments, deleteExperiment, changeExperimentName };