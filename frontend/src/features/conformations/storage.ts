import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {Job, InferenceRequest} from "src/features/conformations/types";
import {
  ConformationsService,
  IntegratorsRequest,
  JobsACommonControllerForJobsManagementService,
  OpenAPI,
  ProteinsService
} from "src/refinedApi/client";
import apiConstants from "src/api/constants";


OpenAPI.BASE = apiConstants.hostname;

const useConformationsStore = defineStore("conformations", {
  actions: {
    async inference(request: InferenceRequest): Promise<{
      job: Job | null,
      errors: string[]
    }> {
      const response = await ConformationsService.runJobApiV1ConformationsJobsRunJobIdPost(
        request.jobId
      )
      const errorResponse = obtainErrorResponse(response);
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

      const protein = await ProteinsService.getProteinApiV1ObjectsProteinsProteinIdGet(response.protein_id);

      return {
        job: {
          id: response.job_id,
          name: response.job_name,
          pdbContent: new File([new Blob([protein!.pdb_content!])], protein!.pdb_name),
          timeline: response.timeline.map(x => {
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

      const protein = await ProteinsService.getProteinApiV1ObjectsProteinsProteinIdGet(
        response.protein_id
      )

      return {
        job: {
          id: response.job_id,
          name: response.job_name,
          pdbContent: new File([new Blob([protein?.pdb_content!])], protein!.pdb_name),
          timeline: response.timeline.map(x => {
            return {
              message: x!.message,
              error: x!.error,
              createdAt: x!.created_at
            };
          }),
          properties: {
            pdbFile: new File([new Blob([protein?.pdb_content!])], protein!.pdb_name),
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
