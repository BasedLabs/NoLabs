import {
  ExperimentsService,
  FoldingService,
  type GetJobMetadataResponse,
  JobsACommonControllerForJobsManagementService,
  type nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse,
  type nolabs__refined__application__use_cases__folding__api_models__JobResponse,
  type ProteinResponse,
  ProteinsService,
  nolabs__refined__application__use_cases__folding__api_models__SetupJobRequest,
  UpdateJobRequest,
  WorkflowService,
  BiobuddyService,
  GenerateMsaService,
  WorkflowSchemaModel_Output,
  AllWorkflowSchemasResponse,
  WorkflowSchemaModel_Input,
  ProteinSearchQuery,
  LigandResponse,
  LigandSearchQuery,
  LigandsService,
  Body_upload_ligand_api_v1_objects_ligands_post,
  BindingPocketsService,
  nolabs__refined__application__use_cases__binding_pockets__api_models__JobResponse,
  nolabs__refined__application__use_cases__binding_pockets__api_models__GetJobStatusResponse,
  DiffdockService,
  nolabs__refined__application__use_cases__diffdock__api_models__JobResponse,
  nolabs__refined__application__use_cases__diffdock__api_models__GetJobStatusResponse,
  nolabs__refined__application__use_cases__diffdock__api_models__SetupJobRequest, OpenAPI,
  GetComponentStateResponse,
  nolabs__refined__application__use_cases__msa_generation__api_models__GetJobStatusResponse,
  nolabs__refined__application__use_cases__msa_generation__api_models__JobResponse,
  Body_upload_protein_api_v1_proteins_post,
  Body_update_protein_api_v1_proteins_patch
} from 'src/refinedApi/client';
import {CancelablePromise, ExperimentMetadataResponse} from "../../refinedApi/client";
import apiConstants from "../../api/constants";

OpenAPI.BASE = apiConstants.hostname;


export function getExperimentsApi(): CancelablePromise<Array<ExperimentMetadataResponse>> {
  return ExperimentsService.experimentsApiV1ExperimentsAllGet();
}

export function createExperimentApi(): CancelablePromise<ExperimentMetadataResponse> {
  return ExperimentsService.createExperimentApiV1ExperimentsPost();
}

export function getProtein(proteinId: string): CancelablePromise<(ProteinResponse | null)> {
  return ProteinsService.getProteinApiV1ProteinsProteinIdGet(proteinId);
}

// Delete an experiment
export function deleteExperimentApi(experimentId: string): CancelablePromise<any> {
  return ExperimentsService.deleteExperimentApiV1ExperimentsExperimentIdDelete(experimentId);
}

export function getFoldingJobApi(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__folding__api_models__JobResponse> {
  return FoldingService.getJobApiV1FoldingJobsJobIdGet(jobId);
}

export function getBindingPocketJobApi(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__binding_pockets__api_models__JobResponse> {
  return BindingPocketsService.getJobApiV1BindingPocketsJobsJobIdGet(jobId);
}

export function getDiffDockJobApi(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__diffdock__api_models__JobResponse> {
  return DiffdockService.getJobApiV1DiffdockJobsJobIdGet(jobId);
}

export function getMsaJobApi(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__msa_generation__api_models__JobResponse> {
  return GenerateMsaService.getJobApiV1MsaGenerationJobsJobIdGet(jobId);
}

export function getFoldingJobStatus(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse> {
  return FoldingService.getJobStatusApiV1FoldingJobsJobIdStatusGet(jobId);
}

export function getBindingPocketJobStatus(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__binding_pockets__api_models__GetJobStatusResponse> {
  return BindingPocketsService.getJobStatusApiV1BindingPocketsJobsJobIdStatusGet(jobId);
}

export function getDiffDockJobStatus(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__diffdock__api_models__GetJobStatusResponse> {
  return DiffdockService.getJobStatusApiV1DiffdockJobsJobIdStatusGet(jobId);
}

export function getMsajobStatus(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__msa_generation__api_models__GetJobStatusResponse> {
  return GenerateMsaService.getJobStatusApiV1MsaGenerationJobsJobIdStatusGet(jobId);
}

export function setupDiffDockJob(job: nolabs__refined__application__use_cases__diffdock__api_models__SetupJobRequest): CancelablePromise<any> {
  return DiffdockService.setupJobApiV1DiffdockJobsPost(job)
}

export function setupFoldingJob(job: nolabs__refined__application__use_cases__folding__api_models__SetupJobRequest): CancelablePromise<nolabs__refined__application__use_cases__folding__api_models__JobResponse> {
  return FoldingService.setupJobApiV1FoldingJobsPost(job);
}

export function startDiffDockJob(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__diffdock__api_models__JobResponse> {
  return DiffdockService.startJobApiV1DiffdockJobsRunJobIdPost(jobId);
}

export function startMsaJob(jobId: string): CancelablePromise<any> {
  return GenerateMsaService.runJobApiV1MsaGenerationJobsRunJobIdPost(jobId);
}

export function changeJobName(jobId: string, newName: string): CancelablePromise<any> {
  const jobRequest = { job_name: newName } as UpdateJobRequest;
  return JobsACommonControllerForJobsManagementService.updateApiV1JobsJobsJobIdPatch(jobId, jobRequest);
}

export function createWorkflow(experimentId: string): CancelablePromise<any> {
  return WorkflowService.createSchemaApiV1WorkflowExperimentIdPost(experimentId);
}

export function getWorkflow(workflowId: string): CancelablePromise<(WorkflowSchemaModel_Output | null)> {
  return WorkflowService.getSchemaApiV1WorkflowWorkflowIdGet(workflowId);
}

export function getExistingWorkflows(experimentId: string): CancelablePromise<AllWorkflowSchemasResponse> {
  return WorkflowService.getSchemaApiV1WorkflowAllExperimentIdGet(experimentId);
}

export function sendWorkflowUpdate(workflow: WorkflowSchemaModel_Input): CancelablePromise<any> {
  return WorkflowService.updateWorkflowSchemaApiV1WorkflowPut(workflow)
}

export function startWorkflow(experimentId: string): CancelablePromise<any> {
  return WorkflowService.startWorkflowApiV1WorkflowWorkflowIdStartPost(experimentId)
}

export function resetWorkflow(workflowId: string): CancelablePromise<any> {
  return WorkflowService.resetWorkflowApiV1WorkflowWorkflowIdResetPost({workflow_id: workflowId});
}

export function getComponentState(componentId: string): CancelablePromise<GetComponentStateResponse> {
  return WorkflowService.getComponentStateApiV1WorkflowComponentComponentIdStateGet(componentId);
}

export function deleteWorkflow(workflowId: string): CancelablePromise<any> {
  return WorkflowService.createSchemaApiV1WorkflowWorkflowIdDelete(workflowId);
}

export function getAllProteins(experimentId: string): CancelablePromise<Array<ProteinResponse>> {
  const searchQuery = {name: '', experiment_id: experimentId} as ProteinSearchQuery;
  return ProteinsService.searchProteinsApiV1ProteinsSearchPost(searchQuery);
}

export function uploadProtein(
  experimentId: string,
  name?: string,
  fasta?: Blob,
  pdb?: Blob
): CancelablePromise<ProteinResponse> {
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

export function updateProteinName(proteinId: string, newName: string): CancelablePromise<ProteinResponse> {
  const newRequest = {
    protein_id: proteinId,
    name: newName
  } as Body_update_protein_api_v1_proteins_patch;
  return ProteinsService.updateProteinApiV1ProteinsPatch(newRequest);
}

export function getAllLigands(experimentId: string): CancelablePromise<Array<LigandResponse>> {
  const searchQuery = {name: '', experiment_id: experimentId} as LigandSearchQuery;
  return LigandsService.searchLigandsApiV1ObjectsLigandsSearchPost(searchQuery);
}

export function uploadLigand(
  experimentId: string,
  name?: string,
  smiles?: Blob,
  sdf?: Blob,
): CancelablePromise<LigandResponse> {
  const uploadLigand = {
    experiment_id: experimentId,
    name: name,
    smiles: smiles,
    sdf: sdf
  } as Body_upload_ligand_api_v1_objects_ligands_post;
  return LigandsService.uploadLigandApiV1ObjectsLigandsPost(uploadLigand);
}

export function deleteLigand(ligandId: string): CancelablePromise<any> {
  return LigandsService.deleteLigandApiV1ObjectsLigandsLigandIdDelete(ligandId);
}

export async function checkBiobuddyEnabled(): Promise<boolean> {
  const response = await BiobuddyService.checkBiobuddyEnabledApiV1BiobuddyCheckBiobuddyEnabledGet();
  return response.enabled;
}

