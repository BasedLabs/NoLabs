import {defineStore} from "pinia";
import {Notify} from "quasar";
import {ErrorCodes} from "src/api/errorTypes";
import {obtainErrorResponse} from "src/api/errorWrapper";
import {
  OpenAPI,
  FoldingService,
  ProteinContentResponse,
  ProteinsService,
  FoldingBackendEnum
} from "src/refinedApi/client";
import apiConstants from "src/api/constants";
import {InferenceRequest, Job} from "src/features/aminoAcid/types";
import {AminoAcid} from "src/features/aminoAcid/folding/types";
import {JobsACommonControllerForJobsManagementService} from "../../../refinedApi/client";

OpenAPI.BASE = apiConstants.hostname;

const useFoldingStore = defineStore("folding", {
  actions: {
    async setupJob(experimentId: string, request: InferenceRequest){
      let proteins: ProteinContentResponse[] = [];
      for(const fasta of request.fastas){
        const protein = await ProteinsService.uploadProteinApiV1ProteinsPost({
          experiment_id: experimentId,
          name: fasta.name,
          fasta: fasta
        });
        proteins.push(protein);
      }

      await FoldingService.setupJobApiV1FoldingJobsPost(
        {
          experiment_id: experimentId,
          job_id: request.jobId,
          job_name: request.jobName,
          protein_ids: proteins.map(p => p.id),
          backend: FoldingBackendEnum.ESMFOLD_LIGHT
        }
      );
    },
    async inference(request: InferenceRequest): Promise<{
      job: Job<AminoAcid> | null,
      errors: string[]
    }> {
      let job = await FoldingService.getJobApiV1FoldingJobsJobIdGet(request.jobId!);

      await this.setupJob(
        job.experiment_id,
        request
      )

      job = await FoldingService.startJobApiV1FoldingJobsRunJobIdPost(job.job_id);
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
        ids: job.protein_ids
      });

      return {
        job: {
          id: job.job_id,
          name: job.job_name,
          aminoAcids: proteins.map(aa => {
            return {
              name: aa.name!,
              sequence: aa.fasta_content!,
              pdbFile: new File([new Blob([aa.pdb_content!])], aa.pdb_name)
            };
          }),
          properties: {
            fastas: request.fastas
          }
        },
        errors: []
      };
    },
    async getJob(jobId: string): Promise<{
      job: Job<AminoAcid> | null,
      errors: string[]
    }> {
      const job = await FoldingService.getJobApiV1FoldingJobsJobIdGet(jobId);
      const errorResponse = obtainErrorResponse(job);
      if (errorResponse) {
        if (errorResponse.error_code === ErrorCodes.job_not_found) {
          return {
            job: {
              id: jobId,
              name: "New job",
              aminoAcids: [],
              properties: {
                fastas: []
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
        ids: job.protein_ids
      })

      return {
        job: {
          id: job.job_id,
          name: job.job_name,
          aminoAcids: proteins.map(aa => {
            return {
              name: aa.name!,
              sequence: aa.fasta_content!,
              pdbFile: new File([new Blob([aa.pdb_content!])], aa.pdb_name)
            };
          }),
          properties: {
            fastas: proteins.map(f =>
              new File([new Blob([f.fasta_content!])], f.fasta_name)
            )
          }
        }, errors: []
      };
    },
    async changeJobName(jobId: string, newName: string) {
      await JobsACommonControllerForJobsManagementService.updateApiV1JobsJobsJobIdPatch(jobId, {
        job_name: newName
      });
    }
  }
});

export default useFoldingStore;
