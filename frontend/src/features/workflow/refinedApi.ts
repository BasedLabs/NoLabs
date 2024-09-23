import {
  ExperimentsService,
  FoldingService,
  JobsACommonControllerForJobsManagementService,
  type ProteinContentResponse,
  ProteinsService,
  UpdateJobRequest,
  WorkflowService,
  BiobuddyService,
  GenerateMsaService,
  AllWorkflowSchemasResponse,
  WorkflowSchema_Input,
  LigandContentResponse,
  LigandsService,
  Body_upload_ligand_api_v1_objects_ligands_post,
  BindingPocketsService,
  DiffdockService,
  nolabs__application__use_cases__diffdock__api_models__JobResponse,
  nolabs__application__use_cases__diffdock__api_models__SetupJobRequest, OpenAPI,
  nolabs__application__use_cases__msa_generation__api_models__JobResponse,
  Body_upload_protein_api_v1_proteins_post,
  Body_update_protein_api_v1_proteins_patch,
  nolabs__application__use_cases__folding__api_models__JobResponse,
  nolabs__application__use_cases__folding__api_models__SetupJobRequest,
  nolabs__application__use_cases__binding_pockets__api_models__JobResponse,
  GetJobMetadataResponse,
  ProteinMetadataResponse,
  LigandMetadataResponse,
  LigandSearchMetadataQuery,
  ProteinSearchMetadataQuery, type WorkflowSchema_Output, type GetComponentResponse, GetJobState
} from 'src/refinedApi/client';
import {CancelablePromise, ExperimentMetadataResponse} from "../../refinedApi/client";
import apiConstants from "../../refinedApi/constants";

OpenAPI.BASE = apiConstants.hostname;


export function getExperimentsApi(): CancelablePromise<Array<ExperimentMetadataResponse>> {
  return ExperimentsService.experimentsApiV1ExperimentsAllGet();
}

export function createExperimentApi(): CancelablePromise<ExperimentMetadataResponse> {
  return ExperimentsService.createExperimentApiV1ExperimentsPost();
}

export function getProteinContent(proteinId: string): CancelablePromise<(ProteinContentResponse | null)> {
  return ProteinsService.getProteinContentApiV1ProteinsProteinIdContentGet(proteinId);
}

export function deleteExperimentApi(experimentId: string): CancelablePromise<any> {
  return ExperimentsService.deleteExperimentApiV1ExperimentsExperimentIdDelete(experimentId);
}

export function getJobMetaData(jobId: string): CancelablePromise<GetJobMetadataResponse> {
  return JobsACommonControllerForJobsManagementService.jobMetadataApiV1JobsJobIdMetadataGet(jobId);
}

export function getFoldingJobApi(jobId: string): CancelablePromise<nolabs__application__use_cases__folding__api_models__JobResponse> {
  return FoldingService.getJobApiV1FoldingJobsJobIdGet(jobId);
}

export function getBindingPocketJobApi(jobId: string): CancelablePromise<nolabs__application__use_cases__binding_pockets__api_models__JobResponse> {
  return BindingPocketsService.getJobApiV1BindingPocketsJobsJobIdGet(jobId);
}

export function getDiffDockJobApi(jobId: string): CancelablePromise<nolabs__application__use_cases__diffdock__api_models__JobResponse> {
  return DiffdockService.getJobApiV1DiffdockJobsJobIdGet(jobId);
}

export function getMsaJobApi(jobId: string): CancelablePromise<nolabs__application__use_cases__msa_generation__api_models__JobResponse> {
  return GenerateMsaService.getJobApiV1MsaGenerationJobsJobIdGet(jobId);
}

export function getJobStatus(jobId: string): CancelablePromise<GetJobState> {
  return WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet(jobId);
}

export function setupDiffDockJob(job: nolabs__application__use_cases__diffdock__api_models__SetupJobRequest): CancelablePromise<any> {
  return DiffdockService.setupJobApiV1DiffdockJobsPost(job)
}

export function setupFoldingJob(job: nolabs__application__use_cases__folding__api_models__SetupJobRequest): CancelablePromise<nolabs__application__use_cases__folding__api_models__JobResponse> {
  return FoldingService.setupJobApiV1FoldingJobsPost(job);
}

export function startDiffDockJob(jobId: string): CancelablePromise<nolabs__application__use_cases__diffdock__api_models__JobResponse> {
  return DiffdockService.startJobApiV1DiffdockJobsRunJobIdPost(jobId);
}

export function startMsaJob(jobId: string): CancelablePromise<any> {
  return GenerateMsaService.runJobApiV1MsaGenerationJobsRunJobIdPost(jobId);
}

export function changeJobName(jobId: string, newName: string): CancelablePromise<any> {
  const jobRequest = { job_name: newName } as UpdateJobRequest;
  return JobsACommonControllerForJobsManagementService.updateApiV1JobsJobIdPatch(jobId, jobRequest);
}

export function createWorkflow(experimentId: string): CancelablePromise<any> {
  return WorkflowService.createWorkflowSchemaApiV1WorkflowExperimentIdPost(experimentId);
}

export function getWorkflow(experimentId: string): CancelablePromise<(WorkflowSchema_Output | null)> {
  return WorkflowService.getSchemaApiV1WorkflowExperimentIdGet(experimentId);
}

export function sendWorkflowUpdate(workflow: WorkflowSchema_Input): CancelablePromise<any> {
  return WorkflowService.updateWorkflowSchemaApiV1WorkflowPut(workflow)
}

export function startWorkflow(experimentId: string): CancelablePromise<any> {
  return WorkflowService.startWorkflowApiV1WorkflowExperimentIdStartPost(experimentId)
}

export function startWorkflowComponent(experimentId:string, componentId: string): CancelablePromise<any> {
  return WorkflowService.startComponentApiV1WorkflowWorkflowIdStartComponentIdPost(experimentId, componentId);
}

export function getComponentState(componentId: string): CancelablePromise<GetComponentResponse> {
  return WorkflowService.getComponentStateApiV1WorkflowComponentComponentIdStateGet(componentId);
}

export function deleteWorkflow(experimentId: string): CancelablePromise<any> {
  return WorkflowService.createSchemaApiV1WorkflowWorkflowIdDelete(experimentId);
}

export function getAllProteinsMetadata(experimentId: string): CancelablePromise<Array<ProteinMetadataResponse>> {
  const searchQuery = {name: '', experiment_id: experimentId} as ProteinSearchMetadataQuery;
  return ProteinsService.searchProteinsApiV1ProteinsSearchMetadataPost(searchQuery);
}

export function uploadProtein(
  experimentId: string,
  name?: string,
  fasta?: Blob,
  pdb?: Blob,
  link?: string
): CancelablePromise<ProteinContentResponse> {
  const uploadProtein = {
    experiment_id: experimentId,
    name: name,
    fasta: fasta,
    pdb: pdb
  } as Body_upload_protein_api_v1_proteins_post;
  return ProteinsService.uploadProteinApiV1ProteinsPost(uploadProtein);
}

export function deleteProtein(proteinId: string): CancelablePromise<any> {
  return ProteinsService.deleteProteinApiV1ProteinsProteinIdDelete(proteinId);
}

export function updateProteinName(proteinId: string, newName: string): CancelablePromise<ProteinContentResponse> {
  const newRequest = {
    protein_id: proteinId,
    name: newName
  } as Body_update_protein_api_v1_proteins_patch;
  return ProteinsService.updateProteinApiV1ProteinsPatch(newRequest);
}

export function getAllLigandsMetadata(experimentId: string): CancelablePromise<Array<LigandMetadataResponse>> {
  const searchQuery = {name: '', experiment_id: experimentId} as LigandSearchMetadataQuery;
  return LigandsService.searchLigandsMetadataApiV1ObjectsLigandsSearchMetadataPost(searchQuery);
}

export function uploadLigand(
  experimentId: string,
  name?: string,
  smiles?: Blob,
  sdf?: Blob,
): CancelablePromise<LigandContentResponse> {
  const uploadLigand = {
    experiment_id: experimentId,
    name: name,
    smiles: smiles,
    sdf: sdf
  } as Body_upload_ligand_api_v1_objects_ligands_post;
  return LigandsService.uploadLigandApiV1ObjectsLigandsPost(uploadLigand);
}


export function getLigandContent(ligandId: string): CancelablePromise<(LigandContentResponse | null)> {
  return LigandsService.getLigandContentApiV1ObjectsLigandsLigandIdContentGet(ligandId);
}

export function deleteLigand(ligandId: string): CancelablePromise<any> {
  return LigandsService.deleteLigandApiV1ObjectsLigandsLigandIdDelete(ligandId);
}

export async function checkBiobuddyEnabled(): Promise<boolean> {
  const response = await BiobuddyService.checkBiobuddyEnabledApiV1BiobuddyCheckBiobuddyEnabledGet();
  return response.enabled;
}

