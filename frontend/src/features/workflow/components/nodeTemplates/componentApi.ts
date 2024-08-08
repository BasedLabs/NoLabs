import {
  BindingPocketsService, BlastService, ConformationsService,
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
    running: false,
    result_valid: true
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
  blast: {
    getJob: BlastService.getJobApiV1BlastJobsJobIdGet,
    executionStatus: BlastService.getJobStatusApiV1BlastJobsJobIdStatusGet
  },
  localisation: {
    getJob: LocalisationService.getJobApiV1LocalisationJobsJobIdGet,
    executionStatus: LocalisationService.getJobStatusApiV1LocalisationJobsJobIdStatusGet
  },
  solubility: {
    getJob: SolubilityService.getJobApiV1SolubilityJobsJobIdGet,
    executionStatus: SolubilityService.getJobStatusApiV1SolubilityJobsJobIdStatusGet
  },
  geneOntology: {
    getJob: GeneOntologyService.getJobApiV1GeneOntologyJobsJobIdGet,
    executionStatus: GeneOntologyService.getJobStatusApiV1GeneOntologyJobsJobIdStatusGet
  },
  conformations: {
    getJob: ConformationsService.getJobApiV1ConformationsJobsJobIdGet,
    executionStatus: ConformationsService.getJobStatusApiV1ConformationsJobsJobIdStatusGet
  },
  proteinDesign: {
    getJob: ProteinDesignService.getJobApiV1ProteinDesignJobsJobIdGet,
    executionStatus: ProteinDesignService.getJobStatusApiV1ProteinDesignJobsJobIdStatusGet
  },
  smallMoleculesDesign: {
    getJob: SmallMoleculesDesignService.getJobApiV1SmallMoleculesDesignJobsJobIdGet,
    executionStatus: SmallMoleculesDesignService.getJobStatusApiV1SmallMoleculesDesignJobsJobIdStatusGet
  }
}

export default componentApi;
