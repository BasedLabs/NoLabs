import {defineStore} from "pinia";
import {OpenAPI, SmallMoleculesDesignService} from "../../api/client";
import apiConstants from "../../api/constants";
import {obtainErrorResponse} from "../../api/errorWrapper";
import {Notify} from "quasar";
import {ExperimentListItem} from "../../components/types";

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
    async getExperiment(experimentId: string): Promise<Experiment> {
      const response = await SmallMoleculesDesignService.getExperimentApiV1SmallMoleculesDesignExperimentExperimentIdGet(
        experimentId
      );

      ensureNotError(response);

      const pdbFile = response.properties.pdb_file ? new File([new Blob([response.properties.pdb_file])], response.properties.pdb_file_name!) : null;

      return {
        id: response.experiment_id,
        name: response.experiment_name,
        samplingAllowed: response.status.sampling_allowed,
        running: response.status.running,
        createdAt: new Date(response.created_at),
        properties: {
          epochs: response.properties.epochs,
          minscore: response.properties.minscore,
          batchSize: response.properties.batch_size,
          centerX: response.properties.center_x,
          centerY: response.properties.center_y,
          centerZ: response.properties.center_z,
          sizeX: response.properties.size_x,
          sizeY: response.properties.size_y,
          sizeZ: response.properties.size_z,
          pdbFile: pdbFile
        }
      };
    },
    async smilesData(experimentId: string): Promise<Smiles[]> {
      const smiles = await SmallMoleculesDesignService.smilesApiV1SmallMoleculesDesignExperimentExperimentIdSmilesGet(
        experimentId
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
    async status(experimentId: string): Promise<Status>{
      const status = await SmallMoleculesDesignService.statusApiV1SmallMoleculesDesignExperimentExperimentIdStatusGet(experimentId);

      ensureNotError(status);

      return {
        running: status.running,
        samplingAllowed: status.sampling_allowed
      }
    },
    async logs(experimentId: string): Promise<Logs> {
      const logs = await SmallMoleculesDesignService.logsApiV1SmallMoleculesDesignExperimentExperimentIdLogsGet(experimentId);

      ensureNotError(logs);

      return {
        output: logs.output,
        dockingOutput: logs.docking_output,
        errors: logs.errors
      }
    },
    async startExperiment(experimentId: string) {
      await SmallMoleculesDesignService.learningApiV1SmallMoleculesDesignExperimentExperimentIdLearningPost(experimentId);
    },
    async startSampling(experimentId: string, numberOfMoleculesToGenerate: number) {
      await SmallMoleculesDesignService.samplingApiV1SmallMoleculesDesignExperimentExperimentIdSamplingPost(experimentId, {
        number: numberOfMoleculesToGenerate
      });
    },
    async stopExperiment(experimentId: string) {
      await SmallMoleculesDesignService.stopApiV1SmallMoleculesDesignExperimentExperimentIdStopPost(experimentId);
    },
    async deleteExperiment(experimentId: string) {
      await SmallMoleculesDesignService.deleteApiV1SmallMoleculesDesignExperimentExperimentIdDelete(experimentId);
    },
    async saveProperties(experimentId: string, properties: Properties): Promise<void> {
      if (!properties.pdbFile) {
        throw new Error('Pdb file is not specified');
      }

      await SmallMoleculesDesignService.savePropertiesApiV1SmallMoleculesDesignExperimentExperimentIdPropsPost(
        experimentId,
        properties.centerX,
        properties.centerY,
        properties.centerZ,
        properties.sizeX,
        properties.sizeY,
        properties.sizeZ,
        {
          pdb_file: properties.pdbFile!
        },
        properties.epochs,
        properties.batchSize,
        properties.minscore
      )
    },
    // Experiments
    async getExperiments(): Promise<{ experiments: ExperimentListItem[] | null, errors: string[] }> {
      const response = await SmallMoleculesDesignService.experimentsApiV1SmallMoleculesDesignExperimentsGet();
      const experiments: ExperimentListItem[] = [];
      for (let i = 0; i < response.length; i++) {
        experiments.push({
          id: response[i].id,
          name: response[i].name
        });
      }
      return {
        experiments: experiments,
        errors: []
      }
    },
    async changeExperimentName(experimentId: string, newName: string) {
      await SmallMoleculesDesignService.changeExperimentNameApiV1SmallMoleculesDesignExperimentNamePost({
        id: experimentId,
        name: newName
      });
    },
    async createExperiment(): Promise<{ experiment: ExperimentListItem | null, errors: string[] }> {
      const response = await SmallMoleculesDesignService.createExperimentApiV1SmallMoleculesDesignExperimentCreatePost();
      return {
        experiment: {
          id: response.id,
          name: response.name
        }, errors: []
      }
    }
  }
});

export default useSmallMoleculesDesignStore;
