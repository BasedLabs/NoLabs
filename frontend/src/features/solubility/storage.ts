import { defineStore } from 'pinia';
import { inference, getExperiment } from './api';
import { AminoAcid } from './types';

export const useSolubilityStore = defineStore('solubility', {
  state: () => ({
    experimentId: null as string | null,
    experimentName: null as string | null,
    aminoAcids: [] as AminoAcid[]
  }),

  actions: {
    async fetchInference(experimentId: string, experimentName: string, aminoAcidSequence: string, fastas: any) {
      try {
        const response = await inference(experimentId, experimentName, aminoAcidSequence, fastas);
        const aminoAcids: AminoAcid[] = [];
        for(const aminoAcid of response.data.amino_acids){
            aminoAcids.push({
              id: aminoAcid.name,
              name: aminoAcid.name,
              sequence: aminoAcid.sequence,
              probability: aminoAcid.soluble_probability
            })
        }
        this.aminoAcids = aminoAcids;
      } catch (error) {
        console.error(error);
        // Handle error
      }
    },

    async fetchExperiment(experimentId: string) {
      try {
        const response = await getExperiment(experimentId);
        // Handle the response, update state, etc.
      } catch (error) {
        console.error(error);
        // Handle error
      }
    },

    // Add other actions as necessary
  },
});