import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {Job, InferenceRequest} from "src/features/conformations/types";
import {
  ConformationsService,
  IntegratorsRequest,
  JobsACommonControllerForJobsManagementService,
  OpenAPI, ProteinsService
} from "src/refinedApi/client";
import apiConstants from "src/api/constants";


OpenAPI.BASE = apiConstants.hostname;

const useConformationsStore = defineStore("conformations", {
  actions: {
    async saveParameters(request: InferenceRequest){
      let job = await ConformationsService.getJobApiV1ConformationsJobsJobIdGet(request.jobId);

      const protein = await ProteinsService.uploadProteinApiV1ProteinsPost({
        experiment_id: job.experiment_id,
        name: request.pdbFile.name,
        pdb: request.pdbFile
      });

      await ConformationsService.setupJobApiV1ConformationsJobsPost(
        {
          protein_id: protein.id,
          job_id: request.jobId,
          job_name: request.jobName,
          total_frames: request.totalFrames,
          temperature_k: request.temperatureK,
          take_frame_every: request.takeFrameEvery,
          step_size: request.stepSize,
          replace_non_standard_residues: request.replaceNonStandardResidues,
          add_missing_atoms: request.addMissingAtoms,
          add_missing_hydrogens: request.addMissingHydrogens,
          friction_coeff: request.frictionCoeff,
          ignore_missing_atoms: request.ignoreMissingAtoms,
          integrator: request.integrator,
          experiment_id: job.experiment_id
        }
      );
    },
    async setupJob(experimentId: string, request: InferenceRequest) {
      const protein = await ProteinsService.uploadProteinApiV1ProteinsPost({
        experiment_id: experimentId,
        name: request.pdbFile.name,
        pdb: request.pdbFile
      });

      await ConformationsService.setupJobApiV1ConformationsJobsPost(
        {
          protein_id: protein.id,
          job_id: request.jobId,
          job_name: request.jobName,
          total_frames: request.totalFrames,
          temperature_k: request.temperatureK,
          take_frame_every: request.takeFrameEvery,
          step_size: request.stepSize,
          replace_non_standard_residues: request.replaceNonStandardResidues,
          add_missing_atoms: request.addMissingAtoms,
          add_missing_hydrogens: request.addMissingHydrogens,
          friction_coeff: request.frictionCoeff,
          ignore_missing_atoms: request.ignoreMissingAtoms,
          integrator: request.integrator,
          experiment_id: experimentId
        }
      );
    },
    async inference(request: InferenceRequest): Promise<{
      job: Job | null,
      errors: string[]
    }> {
      let job = await ConformationsService.getJobApiV1ConformationsJobsJobIdGet(request.jobId);

      await this.setupJob(job.experiment_id, request);

      job = await ConformationsService.runJobApiV1ConformationsJobsRunJobIdPost(
        request.jobId
      )

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

      let protein = null;

      if(job.protein_id){
        protein = await ProteinsService.getProteinApiV1ProteinsProteinIdGet(job.protein_id);
      }

      return {
        job: {
          id: job.job_id,
          name: job.job_name,
          pdbContent: job.md_content ? new File([new Blob([job.md_content!])], protein!.pdb_name) : null,
          timeline: job.timeline.map(x => {
            return {
              message: x.message,
              error: x.error,
              createdAt: x.created_at
            };
          }),
          properties: {
            pdbFile: request.pdbFile,
            totalFrames: request.totalFrames,
            temperatureK: request.temperatureK,
            takeFrameEvery: request.takeFrameEvery,
            stepSize: request.stepSize,
            replaceNonStandardResidues: request.replaceNonStandardResidues,
            addMissingAtoms: request.addMissingAtoms,
            addMissingHydrogens: request.addMissingHydrogens,
            frictionCoeff: request.frictionCoeff,
            ignoreMissingAtoms: request.ignoreMissingAtoms,
            integrator: request.integrator,
          }
        },
        errors: []
      };
    },
    async getJob(jobId: string): Promise<{
      job: Job | null,
      errors: string[]
    }> {
      const response = await ConformationsService.getJobApiV1ConformationsJobsJobIdGet(jobId);
      const errorResponse = obtainErrorResponse(response);
      if (errorResponse) {
        if (errorResponse.error_code === ErrorCodes.job_not_found) {
          return {
            job: {
              id: jobId,
              name: 'New job',
              pdbContent: null,
              timeline: [],
              properties: {
                pdbFile: null,
                totalFrames: 10000,
                temperatureK: 273.15,
                takeFrameEvery: 1000,
                stepSize: 0.002,
                replaceNonStandardResidues: false,
                addMissingAtoms: false,
                addMissingHydrogens: true,
                frictionCoeff: 1.0,
                ignoreMissingAtoms: false,
                integrator: IntegratorsRequest.LANGEVIN_INTEGATOR
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

      let sourceProtein = null;

      if (response.protein_id) {
        sourceProtein = await ProteinsService.getProteinApiV1ProteinsProteinIdGet(
          response.protein_id
        );
      }

      return {
        job: {
          id: response.job_id,
          name: response.job_name,
          pdbContent: response.md_content ? new File([new Blob([response?.md_content!])], sourceProtein!.pdb_name) : null,
          timeline: response.timeline.map(x => {
            return {
              message: x!.message,
              error: x!.error,
              createdAt: x!.created_at
            };
          }),
          properties: {
            pdbFile: sourceProtein ? new File([new Blob([sourceProtein?.pdb_content!])], sourceProtein!.pdb_name) : null,
            totalFrames: response.total_frames,
            temperatureK: response.temperature_k,
            takeFrameEvery: response.take_frame_every,
            stepSize: response.step_size,
            replaceNonStandardResidues: response.replace_non_standard_residues,
            addMissingAtoms: response.add_missing_atoms,
            addMissingHydrogens: response.add_missing_hydrogens,
            frictionCoeff: response.friction_coeff,
            ignoreMissingAtoms: response.ignore_missing_atoms,
            integrator: response.integrator!
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

export default useConformationsStore;
