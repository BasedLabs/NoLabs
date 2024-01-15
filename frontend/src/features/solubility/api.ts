import { AxiosResponse } from 'axios';
import { SolubilityApi, Configuration, RunSolubilityResponse, NolabsApiModelsSolubilityGetExperimentResponse } from '../../api/client';

const configuration = new Configuration(
    {
        basePath: 'http://127.0.0.1:8000'
    }
)

async function inference(experimentId, experimentName, aminoAcidSequence, fastas): Promise<AxiosResponse<RunSolubilityResponse, any>> {
    const api = new SolubilityApi(configuration);
    return await api.inferenceApiV1SolubilityInferencePost(experimentName, experimentId, aminoAcidSequence, fastas);
}

async function getExperiment(experimentId): Promise<AxiosResponse<NolabsApiModelsSolubilityGetExperimentResponse, any>> {
    const api = new SolubilityApi(configuration);
    return await api.getExperimentApiV1SolubilityGetExperimentGet(experimentId);
}

export { inference, getExperiment };