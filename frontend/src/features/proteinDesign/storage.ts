import { defineStore } from 'pinia';
import { changeExperimentName, deleteExperiment, getExperiments, getExperiment, inference } from './api';
import { Experiment, ExperimentListItem, InferenceRequest } from './types';
import { uuidv4 } from '../../utils';

const useProteinDesignStore = defineStore('proteinDesign', {
  state: () => ({
    experiment: null as Experiment,
    experiments: [] as ExperimentListItem[]
  }),
  actions: {
    async inference(request: InferenceRequest) {
      const response = await inference(request.pdbFile, request.contig, request.numberOfDesigns, request.timesteps,
        request.hotspots, request.experimentName, request.experimentId);
      if (response.errors && response.errors.length > 0) {
        throw new Error(response.errors[0]);
      }
      this.experiment = {
        id: response.experiment_id,
        name: response.experiment_name,
        pdbsContent: response.pdbs_content as string[]
      }

      if (!this.experiments.find(x => x?.id === response.experiment_id)) {
        this.experiments.push({
          id: response.experiment_id,
          name: response.experiment_name
        })
      }
    },
    async getExperiment(experimentId: string) {
      const response = await getExperiment(experimentId);
      const experiment = this.experiments.find(x => x?.id === experimentId);
      if ('errors' in response) {
        this.experiment = {
          id: experimentId,
          name: 'New experiment',
          pdbsContent: []
        }
      }
      else {
        this.experiment = {
          id: response.experiment_id,
          name: response.experiment_name,
          pdbsContent: response.pdbs_content
        }
      }
    },
    async getExperiments() {
      const response = await getExperiments();
      const experiments: ExperimentListItem[] = []
      for (let i = 0; i < response.length; i++) {
        experiments.push({
          id: response[i].id,
          name: response[i].name
        })
      }
      this.experiments = experiments;
    },
    async deleteExperiment(experimentId: string) {
      if (this.experiment && this.experiment.id === experimentId) {
        this.experiment = null;
      }
      this.experiments = this.experiments.filter(x => x?.id !== experimentId);
      await deleteExperiment(experimentId);
    },
    async changeExperimentName(experimentId: string, newName: string) {
      const experiment = this.experiments.find(x => x?.id === experimentId);

      if (experiment) {
        experiment.name = newName;
      }
      if (this.experiment && this.experiment.id === experimentId) {
        this.experiment = {
          ...this.experiment,
          name: newName
        }
      }
      await changeExperimentName(experimentId, newName);
    },
    addExperiment() {
      this.experiments.push({
        id: uuidv4(),
        name: 'New experiment'
      });
    }
  }
});

export default useProteinDesignStore;