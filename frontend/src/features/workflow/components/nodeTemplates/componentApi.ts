import {
  DiffdockService,
  FoldingService,
  WorkflowService,
  OpenAPI,
  RfdiffusionService
} from "src/refinedApi/client";
import apiConstants from "src/refinedApi/constants";

OpenAPI.BASE = apiConstants.hostname;

const componentApi = {
  esmfoldLight: {
    getJob: FoldingService.getJobApiV1FoldingJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  esmfold: {
    getJob: FoldingService.getJobApiV1FoldingJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  rosettafold: {
    getJob: FoldingService.getJobApiV1FoldingJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  diffdock: {
    getJob: DiffdockService.getJobApiV1DiffdockJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  proteinDesign: {
    getJob: RfdiffusionService.getJobApiV1RfdiffusionJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  proteinMpnn: {
    getJob: null,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  }
}

export default componentApi;
