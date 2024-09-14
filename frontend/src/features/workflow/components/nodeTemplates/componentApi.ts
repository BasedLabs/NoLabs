import {
  BindingPocketsService, BlastService, ConformationsService,
  DiffdockService,
  FoldingService, GeneOntologyService,
  GenerateMsaService,
  LocalisationService,
  WorkflowService,
  OpenAPI, ProteinDesignService, SmallMoleculesDesignService, SolubilityService
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
  bindingPockets: {
    getJob: BindingPocketsService.getJobApiV1BindingPocketsJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  msaGeneration: {
    getJob: GenerateMsaService.getJobApiV1MsaGenerationJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  diffdock: {
    getJob: DiffdockService.getJobApiV1DiffdockJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  blast: {
    getJob: BlastService.getJobApiV1BlastJobsJobIdGet,
    executionStatus: BlastService.getJobStatusApiV1BlastJobsJobIdStatusGet
  },
  localisation: {
    getJob: LocalisationService.getJobApiV1LocalisationJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  solubility: {
    getJob: SolubilityService.getJobApiV1SolubilityJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  geneOntology: {
    getJob: GeneOntologyService.getJobApiV1GeneOntologyJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  conformations: {
    getJob: ConformationsService.getJobApiV1ConformationsJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  proteinDesign: {
    getJob: ProteinDesignService.getJobApiV1ProteinDesignJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  },
  smallMoleculesDesign: {
    getJob: SmallMoleculesDesignService.getJobApiV1SmallMoleculesDesignJobsJobIdGet,
    executionStatus: WorkflowService.getJobStateApiV1WorkflowJobJobIdStateGet
  }
}

export default componentApi;
