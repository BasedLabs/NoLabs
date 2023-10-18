import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        aminoAcid: {
            prevExperimentName: '',
            currentExperimentName: '',
            experiment: null,
            experiments: []
        },
        drugTarget: {
            prevExperimentName: '',
            currentExperimentName: '',
            experiment: null,
            experiments: []
        }
    },
    actions: {
        aminoAcid: {
            async inference({ commit }, { form }) {
                const response = await axios({
                    method: 'post',
                    url: apiConstants.aminoAcid.inference.path,
                    data: new FormData(form),
                    headers: { "Content-Type": "multipart/form-data" },
                });
                commit(apiConstants.aminoAcid.inference.mutation, response.data);
            },
            async getAllExperiments({ commit }) {
                const response = await axios.get(apiConstants.drugTarget.experiments.path);
                commit(apiConstants.aminoAcid.experiments.mutation, response.data);
            },
            async deleteExperiment({ commit }, { experiment }) {
                return await axios.delete(apiConstants.aminoAcid.deleteExperiment.path, { name: experiment.name });
            },
            async loadExperiment({ commit }, { name }) {
                const response = await axios.get(apiConstants.drugTarget.loadExperiment.path, { params: {name} });
                commit(apiConstants.aminoAcid.loadExperiment.mutation, {name: name, experiment: response.data});
            }
        },
        drugTarget: {
            async drugTargetInference({ commit }, { form }) {
                const response = await axios({
                    method: 'post',
                    url: apiConstants.drugTarget.inference.path,
                    data: new FormData(form),
                    headers: { "Content-Type": "multipart/form-data" },
                });
                commit(apiConstants.drugTarget.inference.mutation, response.data);
            },
            async getAllExperiments({ commit }) {
                const response = await axios.get(apiConstants.drugTarget.experiments.path);
                commit(apiConstants.drugTarget.experiments.mutation, response.data);
            },
            async deleteExperiment({ commit }, { experiment }) {
                return await axios.delete(apiConstants.drugTarget.deleteExperiment.path, { name: experiment.name });
            },
            async loadExperiment({ commit }, { name }) {
                const response = await axios.get(apiConstants.drugTarget.loadExperiment.path, { params: {name} });
                commit(apiConstants.drugTarget.loadExperiment.mutation, {name: name, experiment: response.data});
            }
        }
    },
    mutations: {
        aminoAcid: {
            loadExperiment(state, data) {
                let {name, experiment} = data;
                state.drugTarget.prevExperimentName = state.aminoAcid.currentExperimentName;
    
                if(!experiment || Object.keys(experiment).length == 0){
                    experiment = {name: name}
                }
    
                state.aminoAcid.experiment = experiment;
                state.aminoAcid.currentExperimentName = experiment.name;
            },
            getAllExperiments(state, experiments) {
                if (experiments.length === 0) {
                    api.addExperiment();
                }
                state.aminoAcid.experiments = experiments;
            },
            inference(state, experiment) {
                state.aminoAcid.experiment = experiment;
            },
            addExperiment(state, experiment) {
                state.aminoAcid.experiments.push(experiment);
            }
        },
        drugTarget: {
            loadExperiment(state, data) {
                let {name, experiment} = data;
                state.drugTarget.prevExperimentName = state.drugTarget.currentExperimentName;
    
                if(!experiment || Object.keys(experiment).length == 0){
                    experiment = {name: name}
                }
    
                state.drugTarget.experiment = experiment;
                state.drugTarget.currentExperimentName = experiment.name;
            },
            getAllExperiments(state, experiments) {
                if (experiments.length === 0) {
                    api.addExperiment();
                }
                state.drugTarget.experiments = experiments;
            },
            inference(state, experiment) {
                state.drugTarget.experiment = experiment;
            },
            addExperiment(state, experiment) {
                state.drugTarget.experiments.push(experiment);
            }
        }
    }
});

export const api = {
    aminoAcid: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.aminoAcid.experiments.action);
        },
        inference: async (form) => {
            return await store.dispatch(apiConstants.aminoAcid.inference.action, { form });
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.aminoAcid.deleteExperiment.action, {experiment});
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.aminoAcid.loadExperiment.action, { name: experiment.name });
        },
        addExperiment: () => {
            let latestExperiment = 1;
            for (const experiment of store.state.aminoAcid.experiments) {
                const experimentNumber = parseInt(experiment.name.split(' ')[1]);
                if (experimentNumber >= latestExperiment)
                    latestExperiment = experimentNumber + 1;
            }
            store.dispatch(apiConstants.aminoAcid.addExperiment.action, { experiment: { name: `Experiment ${latestExperiment}` } });
        }
    },
    drugTarget: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.drugTarget.experiments.action);
        },
        inference: async (form) => {
            return await store.dispatch(apiConstants.drugTarget.inference.action, {form});
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.deleteExperiment.action, {experiment});
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.loadExperiment.action, { name: experiment.name });
        },
        addExperiment: () => {
            let latestExperiment = 1;
            for (const experiment of store.state.drugTarget.experiments) {
                const experimentNumber = parseInt(experiment.name.split(' ')[1]);
                if (experimentNumber >= latestExperiment)
                    latestExperiment = experimentNumber + 1;
            }
            store.dispatch(apiConstants.drugTarget.addExperiment.action, { experiment: { name: `Experiment ${latestExperiment}` } });
        }
    }
}

export default store;