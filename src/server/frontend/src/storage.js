import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        aminoAcidData: {
            inference: {}
        },
        drugTargetData: {
            experiment: {},
            experiments: []
        }
    },
    actions: {
        async aminoAcidInference({commit}, { formData }) {
            const response = await axios.post(apiConstants.aminoAcid.inference.path, formData);
            commit(apiConstants.aminoAcid.inference.mutation, response.data);
        },
        async getAllExperiments({commit}) {
            const response = await axios.get(apiConstants.drugTargetDiscovery.experiments.path);
            commit(apiConstants.drugTargetDiscovery.experiments.mutation, response.data);
        },
        async addExperiment({commit}, {experimentId}) {
            return await axios.post(apiConstants.drugTargetDiscovery.addExperiment.path, {id: experimentId});
        },
        async deleteExperiment({commit}, {experimentId}){
            return await axios.delete(apiConstants.drugTargetDiscovery.deleteExperiment.path, {id: experimentId});
        },
        async loadExperiment({commit}, {experimentId}){
            const response = await axios.get(apiConstants.drugTargetDiscovery.loadExperiment.action, {id: experimentId});
            commit(apiConstants.drugTargetDiscovery.loadExperiment.mutation, response.data);
        }
    },
    mutations: {
        aminoAcidInference(state, inference) {
            state.aminoAcidData.inference = inference;
        },
        loadExperiment(state, experiment) {
            state.drugTargetData.experiment = experiment;
        },
        getAllExperiments(state, experiments) {
            state.drugTargetData.experiments = experiments;
        }
    }
});

export const api = {
    aminoAcidLab: {
        inference: async (formData) => {
            return await store.dispatch(apiConstants.aminoAcid.inference.action, {formData});
        }
    },
    drugTargetDiscovery: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.experiments.action);
        },
        addExperiment: async () => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.addExperiment.action);
        },
        deleteExperiment: async () => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.deleteExperiment.action);
        },
        loadExperiment: async (experimentId) => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.loadExperiment.action, {id: experimentId});
        }
    }
}

export default store;