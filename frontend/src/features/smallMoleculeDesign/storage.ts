import {defineStore} from "pinia";
import {
  OpenAPI,
  SmallMoleculesDesignService,
  ProteinsService,
  LigandsService,
  JobsACommonControllerForJobsManagementService, ProteinDesignService
} from "../../refinedApi/client";
import apiConstants from "../../api/constants";
import {obtainErrorResponse} from "../../api/errorWrapper";
import {Notify} from "quasar";

OpenAPI.BASE = apiConstants.hostname;

const ensureNotError = (response: any) => {
  const errorResponse = obtainErrorResponse(response);
  if (errorResponse) {
    for (const error of errorResponse.errors) {
      Notify.create({
        type: "negative",
        closeBtn: 'Close',
        message: error
      });
    }
    throw new Error(errorResponse.errors.join(', '));
  }
}

const useSmallMoleculesDesignStore = defineStore("smallMoleculesDesignStore", {
  actions: {
    async getJob(jobId: string): Promise<Job> {
      const job = await SmallMoleculesDesignService.getJobApiV1SmallMoleculesDesignJobsJobIdGet(
        jobId
      );
      const status = await SmallMoleculesDesignService.getJobStatusApiV1SmallMoleculesDesignJobsJobIdStatusGet(jobId);

      const protein = await ProteinsService.getProteinApiV1ObjectsProteinsProteinIdGet(job!.protein_id);

      ensureNotError(job);
      ensureNotError(status);

      const pdbFile = protein ? new File([new Blob([protein!.pdb_content!])], protein.pdb_name!) : null;

      return {
        id: job!.job_id,
        name: job!.job_name,
        samplingAllowed: status.sampling_allowed,
        running: status.running,
        properties: {
          epochs: job!.epochs,
          minscore: job!.minscore,
          batchSize: job!.batch_size,
          centerX: job!.center_x,
          centerY: job!.center_y,
          centerZ: job!.center_z,
          sizeX: job!.size_x,
          sizeY: job!.size_y,
          sizeZ: job!.size_z,
          pdbFile: pdbFile,
          samplingSize: job!.sampling_size
        }
      };
    },
    async smilesData(jobId: string): Promise<Smiles[]> {
      const smiles = await SmallMoleculesDesignService.getJobSmilesApiV1SmallMoleculesDesignJobsJobIdSmilesGet(
        jobId
      );

      ensureNotError(smiles);

      return smiles.map(x => {
        return {
          smiles: x.smiles,
          drugLikeness: x.drug_likeness,
          score: x.score,
          stage: x.stage,
          createdAt: new Date(x.created_at),
        };
      });
    },
    async status(jobId: string): Promise<Status>{
      const status = await SmallMoleculesDesignService.getJobStatusApiV1SmallMoleculesDesignJobsJobIdStatusGet(jobId);

      ensureNotError(status);

      return {
        running: status.running,
        samplingAllowed: status.sampling_allowed
      }
    },
    async logs(jobId: string): Promise<Logs> {
      const logs = await SmallMoleculesDesignService.getJogLogsApiV1SmallMoleculesDesignJobsJobIdLogsGet(jobId);

      ensureNotError(logs);

      return {
        output: logs.output,
        dockingOutput: logs.docking_output,
        errors: logs.errors
      }
    },
    async startJob(jobId: string) {
      await SmallMoleculesDesignService.runLearningStageJobApiV1SmallMoleculesDesignJobsJobIdRunLearningPost(jobId);
    },
    async startSampling(jobId: string) {
      await SmallMoleculesDesignService.runSamplingStageJobApiV1SmallMoleculesDesignJobsJobIdRunSamplingPost(jobId);
    },
    async stopJob(jobId: string) {
      await SmallMoleculesDesignService.stopJobApiV1SmallMoleculesDesignJobsJobIdStopPost(jobId);
    },
    async saveProperties(jobId: string, properties: Properties): Promise<void> {
      if (!properties.pdbFile) {
        throw new Error('Pdb file is not specified');
      }

      const job = await ProteinDesignService.getJobApiV1ProteinDesignJobsJobIdGet(jobId);

      const protein = await ProteinsService.uploadProteinApiV1ObjectsProteinsPost({
        experiment_id: job.experiment_id,
        name: properties.pdbFile.name,
        pdb: properties.pdbFile
      });

      await SmallMoleculesDesignService.setupJobApiV1SmallMoleculesDesignJobsPost(
        {
          experiment_id: job.experiment_id,
          protein_id: protein.id,
          center_x: properties.centerX,
          center_y: properties.centerY,
          center_z: properties.centerZ,
          size_x: properties.sizeX,
          size_y: properties.sizeY,
          size_z: properties.sizeZ,
          batch_size: properties.batchSize,
          minscore: properties.minscore,
          epochs: properties.epochs,
          job_id: jobId,
          sampling_size: properties.samplingSize
        }
      )
    },
    async changeJobName(jobId: string, newName: string) {
      await JobsACommonControllerForJobsManagementService.updateApiV1JobsJobsJobIdPatch(jobId, {
        job_name: newName
      });
    }
  }
});

export default useSmallMoleculesDesignStore;
