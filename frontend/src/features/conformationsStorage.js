import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        conformations: {
            experiment: {
                metaData: null,
                data: null
            },
            experiments: []
        },
    },
    mutations: {
        conformations_loadExperiment(state, { experiment }) {
            state.conformations.experiment  = experiment;
        },
        conformations_getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            const experimentsArray = [];
            for (const [key, value] of Object.entries(experiments)) {
                experimentsArray.push({
                    metaData: {
                        id: key,
                        name: value.name,
                        progress: value.progress,
                        date: value.date
                    }
                });
            }
            state.conformations.experiments = experimentsArray;
        },
        conformations_inference(state, experiment) {
            state.conformations.experiment = experiment;
        },
        conformations_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.metaData.id = id;
            state.conformations.experiments.push(experiment);
        },
        conformations_deleteExperiment(state, experiment) {
            state.conformations.experiments = state.conformations.experiments.filter(exp => exp.metaData.id !== experiment.metaData.id);
        }
    },
    actions: {
        async conformations_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.generateId.path);
            commit(apiConstants.conformations.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async conformations_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.metaData.id ?? '');
            formData.append('experimentName', experiment.metaData.name);
            const response = await axios({
                method: 'post',
                url: apiConstants.conformations.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.conformations.inference.mutation, response.data);
        },
        async conformations_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.conformations.experiments.path);
            commit(apiConstants.conformations.experiments.mutation, response.data);
        },
        async conformations_deleteExperiment({ commit }, { experiment }) {
            await axios.delete(apiConstants.conformations.deleteExperiment.path, { data: { id: experiment.metaData.id } });
            commit(apiConstants.conformations.deleteExperiment.mutation, experiment);
        },
        async conformations_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.conformations.loadExperiment.path, { params: { id: experiment.metaData.id, name: experiment.metaData.name } });
            commit(apiConstants.conformations.loadExperiment.mutation, { experiment: response.data });
        },
        async conformations_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment)
                await axios.post(apiConstants.conformations.changeExperimentName.path, { id: experiment.id, name: experiment.name });
        }
    },
});

export const api = {
    conformations: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.conformations.experiments.action);
        },
        inference: async (payload) => {
            return await store.dispatch(apiConstants.conformations.inference.action, { payload });
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.conformations.deleteExperiment.action, { experiment });
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.conformations.loadExperiment.action, { experiment });
        },
        addExperiment: () => {
            store.dispatch(apiConstants.conformations.addExperiment.action, { experiment: { metaData: { name: `New experiment` } } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.conformations.changeExperimentName.action, { experiment });
        }
    }
}

export default store;