import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        aminoAcidData: {
            inference: {}
        },
        drugTargetData: {
            prevExperimentName: '',
            currentExperimentName: '',
            experiment: null,
            experiments: []
        }
    },
    actions: {
        addExperiment({ commit }, { experiment }) {
            commit('addExperiment', experiment);
        },
        async drugTargetInference({ commit }, { form }) {
            const response = await axios({
                method: 'post',
                url: apiConstants.drugTargetDiscovery.inference.path,
                data: new FormData(form),
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.drugTargetDiscovery.inference.mutation, response.data);
        },
        async aminoAcidInference({ commit }, { form }) {
            const response = await axios({
                method: 'post',
                url: apiConstants.aminoAcid.inference.path,
                data: new FormData(form),
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.aminoAcid.inference.mutation, response.data);
        },
        async getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.drugTargetDiscovery.experiments.path);
            commit(apiConstants.drugTargetDiscovery.experiments.mutation, response.data);
        },
        async deleteExperiment({ commit }, { experiment }) {
            return await axios.delete(apiConstants.drugTargetDiscovery.deleteExperiment.path, { name: experiment.name });
        },
        async loadExperiment({ commit }, { name }) {
            const response = await axios.get(apiConstants.drugTargetDiscovery.loadExperiment.path, { params: {name} });
            commit(apiConstants.drugTargetDiscovery.loadExperiment.mutation, {name: name, experiment: response.data});
        }
    },
    mutations: {
        aminoAcidInference(state, inference) {
            state.aminoAcidData.inference = inference;
        },
        loadExperiment(state, data) {
            let {name, experiment} = data;
            state.drugTargetData.prevExperimentName = state.drugTargetData.currentExperimentName;

            if(!experiment || Object.keys(experiment).length == 0){
                experiment = {name: name}
            }

            state.drugTargetData.experiment = experiment;
            state.drugTargetData.currentExperimentName = experiment.name;
        },
        getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            state.drugTargetData.experiments = experiments;
        },
        drugTargetInference(state, experiment) {
            state.drugTargetData.experiment = experiment;
        },
        addExperiment(state, experiment) {
            state.drugTargetData.experiments.push(experiment);
        }
    }
});

export const api = {
    aminoAcidLab: {
        inference: async (form) => {
            return await store.dispatch(apiConstants.aminoAcid.inference.action, { form });
        }
    },
    drugTargetDiscovery: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.experiments.action);
        },
        inference: async (form) => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.inference.action, {form});
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.deleteExperiment.action, {experiment});
        },
        loadExperiment: async (experimentName) => {
            return await store.dispatch(apiConstants.drugTargetDiscovery.loadExperiment.action, { name: experimentName });
        },
        addExperiment: () => {
            let latestExperiment = 1;
            for (const experiment of store.state.drugTargetData.experiments) {
                const experimentNumber = parseInt(experiment.name.split(' ')[1]);
                if (experimentNumber >= latestExperiment)
                    latestExperiment = experimentNumber + 1;
            }
            store.dispatch('addExperiment', { experiment: { name: `Experiment ${latestExperiment}` } });
        }
    }
}

export default store;