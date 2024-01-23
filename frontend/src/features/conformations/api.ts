import {
  CancelablePromise,
  nolabs__api_models__protein_design__ExperimentMetadataResponse,
  nolabs__api_models__protein_design__GetExperimentResponse,
  OpenAPI,
  ProteinDesignService,
  RunProteinDesignResponse
} from "src/api/client";
import apiConstants from "src/api/constants";


OpenAPI.BASE = apiConstants.hostname;

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

function createExperiment(): CancelablePromise<nolabs__api_models__protein_design__ExperimentMetadataResponse> {
  return ProteinDesignService.createExperimentApiV1ProteinDesignCreateExperimentGet();
}

function getExperiment(experimentId: string): CancelablePromise<nolabs__api_models__protein_design__GetExperimentResponse> {
  return ProteinDesignService.getExperimentApiV1ProteinDesignExperimentGet(experimentId);
}

function getExperiments(): CancelablePromise<Array<nolabs__api_models__protein_design__ExperimentMetadataResponse>> {
  return ProteinDesignService.experimentsApiV1ProteinDesignExperimentsMetadataGet();
}

function deleteExperiment(experimentId: string): CancelablePromise<any> {
  return ProteinDesignService.deleteExperimentApiV1ProteinDesignExperimentDelete(experimentId);
}

function changeExperimentName(experimentId: string, experimentName: string): CancelablePromise<any> {
  return ProteinDesignService.changeExperimentNameApiV1ProteinDesignChangeExperimentNamePost({
    id: experimentId,
    name: experimentName
  })
}

export { inference, getExperiment, getExperiments, deleteExperiment, changeExperimentName, createExperiment };
