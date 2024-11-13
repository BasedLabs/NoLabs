import {defineStore} from "pinia";
import {Job, InferenceRequest} from "src/features/proteinDesign/types";
import {Notify} from "quasar";
import {ErrorCodes} from "src/refinedApi/errorTypes";
import {obtainErrorResponse} from "src/refinedApi/errorWrapper";
import {
  JobsACommonControllerForJobsManagementService,
  OpenAPI,
  ProteinDesignService,
  ProteinsService
} from "src/refinedApi/client";
import apiConstants from "../../refinedApi/constants";

OpenAPI.BASE = apiConstants.hostname;

const useProteinDesignStore = defineStore("proteinDesign", {
  actions: {
    async setupJob(experimentId: string, request: InferenceRequest) {
      const protein = await ProteinsService.uploadProteinApiV1ProteinsPost({
        experiment_id: experimentId,
        name: request.pdbFile.name,
        pdb: request.pdbFile
      });

      await ProteinDesignService.setupJobApiV1ProteinDesignJobsPost(
        {
          experiment_id: experimentId,
          protein_id: protein.id,
          contig: request.contig,
          number_of_designs: request.numberOfDesigns,
          timesteps: request.timesteps,
          hotspots: request.hotspots,
          job_id: request.jobId,
          job_name: request.jobName
        }
      );
    },
    async inference(request: InferenceRequest): Promise<{
      job: Job | null,
      errors: string[]
    }> {
      let job = await ProteinDesignService.getJobApiV1ProteinDesignJobsJobIdGet(request.jobId!);

      await this.setupJob(job.experiment_id, request);

      job = await ProteinDesignService.runJobApiV1ProteinDesignJobsRunJobIdPost(
        request.jobId!
      );
      const errorResponse = obtainErrorResponse(job);
      if (errorResponse) {
        for (const error of errorResponse.errors) {
          Notify.create({
            type: "negative",
            closeBtn: 'Close',
            message: error
          });
        }
        return {job: null, errors: errorResponse.errors};
      }
      const proteins = await ProteinsService.searchProteinsApiV1ProteinsSearchContentPost({
        ids: job.binder_ids
      })
      return {
        job: {
          id: job.job_id,
          experimentId: job.experiment_id,
          name: job.job_name,
          generatedPdbs: proteins.map(protein => new File([new Blob([protein.pdb_content!])], protein.pdb_name)),
          properties: {
            inputPdbFile: request.pdbFile,
            contig: request.contig,
            numberOfDesigns: request.numberOfDesigns,
            timesteps: request.timesteps,
            hotspots: request.hotspots
          }
        },
        errors: []
      };
    },
    async getJob(jobId: string): Promise<{
      job: Job | null,
      errors: string[]
    }> {
      const response = await ProteinDesignService.getJobApiV1ProteinDesignJobsJobIdGet(jobId);
      const errorResponse = obtainErrorResponse(response);
      if (errorResponse) {
        if (errorResponse.error_code === ErrorCodes.job_not_found) {
          return {
            job: {
              id: jobId,
              experimentId: response.experiment_id,
              name: "New job",
              generatedPdbs: [],
              properties: {
                inputPdbFile: null,
                contig: '',
                numberOfDesigns: 2,
                timesteps: 50,
                hotspots: ''
              }
            }, errors: []
          };
        } else {
          Notify.create({
            type: "negative",
            message: errorResponse.errors[0]
          });
        }

        return {job: null, errors: errorResponse.errors};
      }

      const proteins = await ProteinsService.searchProteinsApiV1ProteinsSearchContentPost({
        ids: response.binder_ids
      });

      const protein = await ProteinsService.getProteinContentApiV1ProteinsProteinIdContentGet(
        response.protein_id
      )

      return {
        job: {
          id: response.job_id,
          experimentId: response.experiment_id,
          name: response.job_name,
          generatedPdbs: proteins.map(p => new File([new Blob([p.pdb_content!])], p.pdb_name)),
          properties: {
            inputPdbFile: new File([new Blob([protein?.pdb_content!], {
              type: 'text/plain'
            })], protein?.pdb_name!),
            contig: response.contig,
            numberOfDesigns: response.number_of_designs!,
            timesteps: response.timesteps ?? 50,
            hotspots: response.hotspots ?? ''
          }
        }, errors: []
      };
    },
    async changeJobName(jobId: string, newName: string) {
      await JobsACommonControllerForJobsManagementService.updateApiV1JobsJobsJobIdPatch(jobId, {
        job_name: newName
      });
    },
  }
});

export default useProteinDesignStore;
