import {defineStore} from 'pinia';

import {
  ExperimentMetadataResponse, OpenAPI,
} from 'src/refinedApi/client';
import {ExperimentListItem} from "src/components/types";
import {ExperimentsService} from "src/refinedApi/client";
import apiConstants from "../refinedApi/constants";

OpenAPI.BASE = apiConstants.hostname;

export const useExperimentsStore = defineStore('experiments', {
  state: () => ({
    experiments: [] as ExperimentMetadataResponse[],
    currentExperiment: null as ExperimentMetadataResponse | null
  }),
  actions: {
    async getExperiments(): Promise<ExperimentListItem[]> {
      const response = await ExperimentsService.experimentsApiV1ExperimentsAllGet();
      this.experiments = response;
      const experiments: ExperimentListItem[] = [];
      for (let i = 0; i < response.length; i++) {
        experiments.push({
          id: response[i].id,
          name: response[i].name
        });
      }
      return experiments;
    },
    async createExperiment(): Promise<ExperimentListItem> {
      const response = await ExperimentsService.createExperimentApiV1ExperimentsPost();
      this.experiments.push(response);
      return {
        id: response.id,
        name: response.name
      };
    },
    async deleteExperiment(experimentId: string) {
      try {
        await ExperimentsService.deleteExperimentApiV1ExperimentsExperimentIdDelete(experimentId);
        this.experiments = this.experiments.filter(exp => exp.id !== experimentId);
      } catch (error) {
        console.error('Error deleting experiment:', error);
      }
    },
    async changeExperimentName(experimentId: string, experimentName: string) {
      try {
        return await ExperimentsService.updateExperimentApiV1ExperimentsPatch({id: experimentId, name: experimentName});
      } catch (error) {
        console.error('Error deleting experiment:', error);
      }
    }
  }
});

export default useExperimentsStore;
