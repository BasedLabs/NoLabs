import { defineStore } from 'pinia';
import { inference, getExperiment } from './api';
import { AminoAcid } from './types';

export const useSolubilityStore = defineStore('solubility', {
  state: () => ({
    experimentId: null as string | null,
    experimentName: null as string | null,
    aminoAcids: [] as AminoAcid[]
  })
});