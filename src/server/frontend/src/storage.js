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
        aminoAcid_addExperiment({ commit }, { experiment }) {
            commit(apiConstants.aminoAcid.addExperiment.mutation, experiment);
        },
        async aminoAcid_inference({ commit }, { payload }) {
            const {form, experiment} = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.id);
            const response = await axios({
                method: 'post',
                url: apiConstants.aminoAcid.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.aminoAcid.inference.mutation, response.data);
        },
        async aminoAcid_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.drugTarget.experiments.path);
            commit(apiConstants.aminoAcid.experiments.mutation, response.data);
        },
        async aminoAcid_deleteExperiment({ commit }, { experiment }) {
            return await axios.delete(apiConstants.aminoAcid.deleteExperiment.path, { id: experiment.id });
        },
        async aminoAcid_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperiment.path, { params: { id: experiment.id } });
            commit(apiConstants.aminoAcid.loadExperiment.mutation, { name: experiment.name, experiment: response.data });
        },
        drugTarget_addExperiment({ commit }, { experiment }) {
            commit(apiConstants.drugTarget.addExperiment.mutation, experiment);
        },
        async drugTarget_inference({ commit }, { payload }) {
            const {form, experiment} = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.id);
            const response = await axios({
                method: 'post',
                url: apiConstants.drugTarget.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.drugTarget.inference.mutation, response.data);
        },
        async drugTarget_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.drugTarget.experiments.path);
            commit(apiConstants.drugTarget.experiments.mutation, response.data);
        },
        async drugTarget_deleteExperiment({ commit }, { experiment }) {
            return await axios.delete(apiConstants.drugTarget.deleteExperiment.path, { id: experiment.id });
        },
        async drugTarget_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.drugTarget.loadExperiment.path, { params: { id: experiment.id } });
            commit(apiConstants.drugTarget.loadExperiment.mutation, { name: experiment.name, experiment: response.data });
        }
    },
    mutations: {
        aminoAcid_loadExperiment(state, data) {
            let { name, experiment } = data;
            state.aminoAcid.prevExperimentName = state.aminoAcid.currentExperimentName;

            if (!experiment || Object.keys(experiment).length == 0) {
                experiment = { name: name }
            }

            state.aminoAcid.experiment = experiment;
            state.aminoAcid.currentExperimentName = experiment.name;
        },
        aminoAcid_getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            state.aminoAcid.experiments = experiments;
        },
        aminoAcid_inference(state, experiment) {
            state.aminoAcid.experiment = experiment;
        },
        aminoAcid_addExperiment(state, experiment) {
            state.aminoAcid.experiments.push(experiment);
        },
        drugTarget_loadExperiment(state, data) {
            let { name, experiment } = data;
            state.drugTarget.prevExperimentName = state.drugTarget.currentExperimentName;

            if (!experiment || Object.keys(experiment).length == 0) {
                experiment = { name: name }
            }

            state.drugTarget.experiment = experiment;
            state.drugTarget.currentExperimentName = experiment.name;
        },
        drugTarget_getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            state.drugTarget.experiments = experiments;
        },
        drugTarget_inference(state, experiment) {
            state.drugTarget.experiment = experiment;
        },
        drugTarget_addExperiment(state, experiment) {
            state.drugTarget.experiments.push(experiment);
        }
    }
});

export const api = {
    aminoAcid: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.aminoAcid.experiments.action);
        },
        inference: async (payload) => {
            return await store.dispatch(apiConstants.aminoAcid.inference.action, { payload });
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.aminoAcid.deleteExperiment.action, { experiment });
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.aminoAcid.loadExperiment.action, { experiment });
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
        inference: async (payload) => {
            return await store.dispatch(apiConstants.drugTarget.inference.action, payload);
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.deleteExperiment.action, { experiment });
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.loadExperiment.action, { experiment });
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