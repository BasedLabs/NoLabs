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
  WorkflowSchemaModel_Output,
  AllWorkflowSchemasResponse
} from 'src/refinedApi/client';
import {CancelablePromise, ExperimentMetadataResponse} from "../../refinedApi/client";


export function getExperimentsApi(): CancelablePromise<Array<ExperimentMetadataResponse>> {
  return ExperimentsService.experimentsApiV1ExperimentsExperimentsAllGet();
}

export function createExperimentApi(): CancelablePromise<ExperimentMetadataResponse> {
  return ExperimentsService.createExperimentApiV1ExperimentsExperimentsPost();
}

// Delete an experiment
export function deleteExperimentApi(experimentId: string): CancelablePromise<any> {
  return ExperimentsService.deleteExperimentApiV1ExperimentsExperimentsExperimentIdDelete(experimentId);
}

export function getFoldingJobsApi(experimentId: string): CancelablePromise<Array<GetJobMetadataResponse>> {
  return JobsACommonControllerForJobsManagementService.jobsMetadataApiV1JobsJobsMetadataGet(experimentId);
}

export function getFoldingJobApi(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__folding__api_models__JobResponse> {
  return FoldingService.getJobApiV1FoldingJobsJobIdGet(jobId);
}

export function getProtein(proteinId: string): CancelablePromise<(ProteinResponse | null)> {
  return ProteinsService.getProteinApiV1ObjectsProteinsProteinIdGet(proteinId);
}

export function getFoldingJobStatus(jobId: string): CancelablePromise<nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse> {
  return FoldingService.getJobStatusApiV1FoldingJobsJobIdStatusGet(jobId);
}

export function setupFoldingJob(job: nolabs__refined__application__use_cases__folding__api_models__SetupJobRequest): CancelablePromise<nolabs__refined__application__use_cases__folding__api_models__JobResponse> {
  return FoldingService.setupJobApiV1FoldingJobsPost(job);
}

export function changeJobName(jobId: string, newName: string): CancelablePromise<any> {
  const jobRequest = { job_name: newName } as UpdateJobRequest;
  return JobsACommonControllerForJobsManagementService.updateApiV1JobsJobsJobIdPatch(jobId, jobRequest);
}

export function getWorkflow(workflowId: string): CancelablePromise<(WorkflowSchemaModel_Output | null)> {
  return WorkflowService.getSchemaApiV1WorkflowWorkflowIdGet(workflowId);
}

export function getExistingWorkflows(experimentId: string): CancelablePromise<AllWorkflowSchemasResponse> {
  return WorkflowService.getSchemaApiV1WorkflowAllExperimentIdGet(experimentId);
}

