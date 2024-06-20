import {
  BindingPocketsService, ConformationsService,
  DiffdockService,
  FoldingService, GeneOntologyService,
  GenerateMsaService,
  LocalisationService,
  OpenAPI, ProteinDesignService, SmallMoleculesDesignService, SolubilityService
} from "src/refinedApi/client";
import apiConstants from "src/refinedApi/constants";

OpenAPI.BASE = apiConstants.hostname;

const mockExecutionStatusApi = () => {
  // TODO remove this mock
  return {
    running: false
  }
}

const componentApi = {
  esmfoldLight: {
    getJob: FoldingService.getJobApiV1FoldingJobsJobIdGet,
    executionStatus: FoldingService.getJobStatusApiV1FoldingJobsJobIdStatusGet
  },
  esmfold: {
    getJob: FoldingService.getJobApiV1FoldingJobsJobIdGet,
    executionStatus: FoldingService.getJobStatusApiV1FoldingJobsJobIdStatusGet
  },
  rosettafold: {
    getJob: FoldingService.getJobApiV1FoldingJobsJobIdGet,
    executionStatus: FoldingService.getJobStatusApiV1FoldingJobsJobIdStatusGet
  },
  bindingPockets: {
    getJob: BindingPocketsService.getJobApiV1BindingPocketsJobsJobIdGet,
    executionStatus: BindingPocketsService.getJobStatusApiV1BindingPocketsJobsJobIdStatusGet
  },
  msaGeneration: {
    getJob: GenerateMsaService.getJobApiV1MsaGenerationJobsJobIdGet,
    executionStatus: GenerateMsaService.getJobStatusApiV1MsaGenerationJobsJobIdStatusGet
  },
  diffdock: {
    getJob: DiffdockService.getJobApiV1DiffdockJobsJobIdGet,
    executionStatus: DiffdockService.getJobStatusApiV1DiffdockJobsJobIdStatusGet
  },
  localisation: {
    getJob: LocalisationService.getJobApiV1LocalisationJobsJobIdGet,
    executionStatus: mockExecutionStatusApi
  },
  solubility: {
    getJob: SolubilityService.getJobApiV1SolubilityJobsJobIdGet,
    executionStatus: mockExecutionStatusApi
  },
  geneOntology: {
    getJob: GeneOntologyService.getJobApiV1GeneOntologyJobsJobIdGet,
    executionStatus: mockExecutionStatusApi
  },
  conformations: {
    getJob: ConformationsService.getJobApiV1ConformationsJobsJobIdGet,
    executionStatus: mockExecutionStatusApi
  },
  proteinDesign: {
    getJob: ProteinDesignService.getJobApiV1ProteinDesignJobsJobIdGet,
    executionStatus: mockExecutionStatusApi
  },
  smallMoleculesDesign: {
    getJob: SmallMoleculesDesignService.getJobApiV1SmallMoleculesDesignJobsJobIdGet,
    executionStatus: SmallMoleculesDesignService.getJobStatusApiV1SmallMoleculesDesignJobsJobIdStatusGet
  }
}

export default componentApi;
