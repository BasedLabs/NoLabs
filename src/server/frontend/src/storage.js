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
        addExperiment({commit}, {experiment}){
            commit('addExperiment', experiment);
        },
        async aminoAcidInference({commit}, { formData }) {
            const response = await axios.post(apiConstants.aminoAcid.inference.path, formData);
            commit(apiConstants.aminoAcid.inference.mutation, response.data);
        },
        async getAllExperiments({commit}) {
            const response = await axios.get(apiConstants.drugTargetDiscovery.experiments.path);
            commit(apiConstants.drugTargetDiscovery.experiments.mutation, response.data);
        },
        async drugTargetInference({commit}, {name}) {
            const response = await axios.post(apiConstants.drugTargetDiscovery.inference.path, {name: name});
            commit(apiConstants.drugTargetDiscovery.loadExperiment.mutation, response.data);
        },
        async deleteExperiment({commit}, {name}){
            return await axios.delete(apiConstants.drugTargetDiscovery.deleteExperiment.path, {name: name});
        },
        async loadExperiment({commit}, {name}){
            const response = await axios.get(apiConstants.drugTargetDiscovery.loadExperiment.action, {name});
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
            if(experiments.length === 0){
                api.addExperiment();
            }
            state.drugTargetData.experiments = experiments;
        },
        drugTargetInference(state, experiment) {
            state.drugTargetData.experiment = experiment;
        },
        addExperiment(state, experiment){
            console.log(state.drugTargetData.experiments.length);
            state.drugTargetData.experiment = experiment;
            state.drugTargetData.experiments.push(experiment);
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
        inference: async () => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.inference.action);
        },
        deleteExperiment: async () => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.deleteExperiment.action);
        },
        loadExperiment: async (experimentId) => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.loadExperiment.action, {id: experimentId});
        },
        addExperiment: () => {
            let latestExperiment = 1;
            for(const experiment of store.state.drugTargetData.experiments){
                const experimentNumber = parseInt(experiment.name.split(' ')[1]);
                if(experimentNumber > latestExperiment)
                    latestExperiment = experimentNumber
            }
            store.dispatch('addExperiment', {name: `Experiment ${latestExperiment}`});
        }
    }
}

export default store;