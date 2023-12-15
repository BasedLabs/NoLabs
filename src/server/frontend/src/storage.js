import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        aminoAcid: {
            experiment: {'metaData': null,
                         'data': null},
            experiments: []
        },
        drugTarget: {
            experiment: null,
            experiments: []
        },
        conformations: {
            experiment: null,
            experiments: []
        }
    },
    mutations: {
        aminoAcid_loadExperiment(state, { experiment }) {
            state.aminoAcid.experiment.metaData = experiment;
        },
        aminoAcid_loadResults(state, { experimentData }) {
            debugger;
            state.aminoAcid.experiment.data = experimentData.data;
        },
        aminoAcid_loadExperimentProgress(state, { experimentProgress }) {
            state.aminoAcid.experiment.metaData.progress = experimentProgress.data;
        },
        aminoAcid_loadExperimentInstanceProgress(state, { experimentInstanceProgress }) {
            state.aminoAcid.experiment.metaData.proteinIds[experimentInstanceProgress.id] = experimentInstanceProgress.progress;
        },
        aminoAcid_getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            const experimentsArray = [];
            for (const [key, value] of Object.entries(experiments)) {
                experimentsArray.push({
                    metaData: {id: key, 
                        name: value.name, 
                        progress: value.progress,
                        date: value.date 
                    }
                });
            }
            state.aminoAcid.experiments = experimentsArray;
        },
        aminoAcid_inference(state, experiment) {
            state.aminoAcid.experiment = experiment;
        },
        aminoAcid_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.metaData.id = id;
            state.aminoAcid.experiments.push(experiment);
        },
        aminoAcid_deleteExperiment(state, experiment) {
            if (state.aminoAcid.experiment && state.aminoAcid.experiment.metaData.id === experiment.metaData.id) {
                state.aminoAcid.experiment = null;
            }
            state.aminoAcid.experiments = state.aminoAcid.experiments.filter(exp => exp.id !== experiment.metaData.id);
        },
        drugTarget_loadExperiment(state, { experiment }) {
            state.drugTarget.experiment = experiment;
        },
        drugTarget_getAllExperiments(state, experiments) {
            const experimentsArray = [];
            for (const [key, value] of Object.entries(experiments)) {
                experimentsArray.push({ id: key, name: value });
            }
            state.drugTarget.experiments = experimentsArray;
        },
        drugTarget_inference(state, experiment) {
            state.drugTarget.experiment = experiment;
        },
        drugTarget_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.id = id;
            state.drugTarget.experiments.push(experiment);
        },
        drugTarget_deleteExperiment(state, experiment) {
            if (state.drugTarget.experiment && state.drugTarget.experiment.id === experiment.id) {
                state.drugTarget.experiment = null;
            }
            state.drugTarget.experiments = state.drugTarget.experiments.filter(exp => exp.id !== experiment.id);
        },
        conformations_loadExperiment(state, { experiment }) {
            state.conformations.experiment = experiment;
        },
        conformations_getAllExperiments(state, experiments) {
            const experimentsArray = [];
            for (const [key, value] of Object.entries(experiments)) {
                experimentsArray.push({ id: key, name: value });
            }
            state.conformations.experiments = experimentsArray;
        },
        conformations_inference(state, experiment) {
            state.conformations.experiment = experiment;
        },
        conformations_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.id = id;
            state.conformations.experiments.push(experiment);
        },
        conformations_deleteExperiment(state, experiment) {
            if (state.conformations.experiment && state.conformations.experiment.id === experiment.id) {
                state.conformations.experiment = null;
            }
            state.conformations.experiments = state.conformations.experiments.filter(exp => exp.id !== experiment.id);
        },
    },
    actions: {
        async aminoAcid_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.generateId.path);
            commit(apiConstants.aminoAcid.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async aminoAcid_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.metaData.id ?? '');
            formData.append('experimentName', experiment.metaData.name);
            const response = await axios({
                method: 'post',
                url: apiConstants.aminoAcid.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
        },
        async aminoAcid_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.aminoAcid.experiments.path);
            commit(apiConstants.aminoAcid.experiments.mutation, response.data);
        },
        async aminoAcid_deleteExperiment({ commit }, { experiment }) {
            await axios.delete(apiConstants.aminoAcid.deleteExperiment.path, { data: { id: experiment.metaData.id } });
            commit(apiConstants.aminoAcid.deleteExperiment.mutation, experiment);
        },
        async aminoAcid_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperiment.path, { params: { id: experiment.metaData.id, name: experiment.metaData.name } });
            commit(apiConstants.aminoAcid.loadExperiment.mutation, { experiment: response.data });
        },
        async aminoAcid_loadResults({ commit }, { experiment, proteinId }) {
            const response = await axios.get(apiConstants.aminoAcid.loadResults.path, 
                { params: 
                    { id: experiment.metaData.id,
                      name: experiment.metaData.name,
                      proteinId: proteinId } });
            commit(apiConstants.aminoAcid.loadResults.mutation, { experimentData: response.data });
        },
        async aminoAcid_loadExperimentProgress({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperimentProgress.path, 
                { params: { id: experiment.metaData.id } });
            commit(apiConstants.aminoAcid.loadExperimentProgress.mutation, { experimentData: response.data });
        },
        async aminoAcid_loadExperimentInstanceProgress({ commit }, { experiment, proteinId }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperimentInstanceProgress.path, 
                { params: 
                    { id: experiment.metaData.id,
                      proteinId: proteinId } });
            commit(apiConstants.aminoAcid.loadExperimentInstanceProgress.mutation, { experimentInstanceProgress: response.data });
        },
        async aminoAcid_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment.metaData)
                await axios.post(apiConstants.aminoAcid.changeExperimentName.path, { id: experiment.metaData.id, name: experiment.metaData.name });
        },
        async drugTarget_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.generateId.path);
            commit(apiConstants.drugTarget.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async drugTarget_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.id ?? '');
            formData.append('experimentName', experiment.name);
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
            await axios.delete(apiConstants.drugTarget.deleteExperiment.path, { data: { id: experiment.id } });
            commit(apiConstants.drugTarget.deleteExperiment.mutation, experiment);
        },
        async drugTarget_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.drugTarget.loadExperiment.path, { params: { id: experiment.id, name: experiment.name } });
            commit(apiConstants.drugTarget.loadExperiment.mutation, { experiment: response.data });
        },
        async drugTarget_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment)
                await axios.post(apiConstants.drugTarget.changeExperimentName.path, { id: experiment.id, name: experiment.name });
        },
        async conformations_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.generateId.path);
            commit(apiConstants.conformations.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async conformations_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.id ?? '');
            formData.append('experimentName', experiment.name);
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
            await axios.delete(apiConstants.conformations.deleteExperiment.path, { data: { id: experiment.id } });
            commit(apiConstants.conformations.deleteExperiment.mutation, experiment);
        },
        async conformations_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.conformations.loadExperiment.path, { params: { id: experiment.id, name: experiment.name } });
            commit(apiConstants.conformations.loadExperiment.mutation, { experiment: response.data });
        },
        async conformations_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment)
                await axios.post(apiConstants.conformations.changeExperimentName.path, { id: experiment.id, name: experiment.name });
        }
    },
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
        loadResults: async (experiment, proteinId) => {
            return await store.dispatch(apiConstants.aminoAcid.loadResults.action, { experiment, proteinId });
        },
        loadExperimentProgress: async (experiment) => {
            return await store.dispatch(apiConstants.aminoAcid.loadExperimentProgress.action, { experiment });
        },
        loadExperimentInstanceProgress: async (experiment, proteinId) => {
            return await store.dispatch(apiConstants.aminoAcid.loadExperimentInstanceProgress.action, { experiment, proteinId });
        },
        addExperiment: () => {
            store.dispatch(apiConstants.aminoAcid.addExperiment.action, { experiment: { name: `New experiment` } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.aminoAcid.changeExperimentName.action, { experiment });
        }
    },
    drugTarget: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.drugTarget.experiments.action);
        },
        inference: async (payload) => {
            return await store.dispatch(apiConstants.drugTarget.inference.action, { payload });
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.deleteExperiment.action, { experiment });
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.loadExperiment.action, { experiment });
        },
        addExperiment: () => {
            store.dispatch(apiConstants.drugTarget.addExperiment.action, { experiment: { name: 'New experiment' } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.drugTarget.changeExperimentName.action, { experiment });
        },
        downloadCombinedPdb: async (experiment, selectedIndex) => {
            return await axios.post(apiConstants.drugTarget.downloadCombinedPdb.path,
                { experiment_id: experiment.id, selected_index: selectedIndex });
        }
    },
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
            store.dispatch(apiConstants.conformations.addExperiment.action, { experiment: { name: 'New experiment' } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.conformations.changeExperimentName.action, { experiment });
        }
    }
}

export default store;