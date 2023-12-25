import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        aminoAcid: {
            experiment: null,
            experiments: []
        },
        drugTarget: {
            experiment: null,
            experiments: []
        },
        conformations: {
            experiment: null,
            experiments: []
        },
        proteinDesign: {
            experiment: null,
            experiments: []
        }
    },
    mutations: {
        aminoAcid_loadExperiment(state, { experiment }) {
            state.aminoAcid.experiment = experiment;
        },
        aminoAcid_getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            const experimentsArray = [];
            for (const [key, value] of Object.entries(experiments)) {
                experimentsArray.push({ id: key, name: value });
            }
            state.aminoAcid.experiments = experimentsArray;
        },
        aminoAcid_inference(state, experiment) {
            state.aminoAcid.experiment = experiment;
        },
        aminoAcid_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.id = id;
            state.aminoAcid.experiments.push(experiment);
        },
        aminoAcid_deleteExperiment(state, experiment) {
            if (state.aminoAcid.experiment && state.aminoAcid.experiment.id === experiment.id) {
                state.aminoAcid.experiment = null;
            }
            state.aminoAcid.experiments = state.aminoAcid.experiments.filter(exp => exp.id !== experiment.id);
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
        proteinDesign_loadExperiment(state, { experiment }) {
            state.proteinDesign.experiment = experiment;
        },
        proteinDesign_getAllExperiments(state, experiments) {
            if (experiments.length === 0) {
                api.addExperiment();
            }
            const experimentsArray = [];
            for (const [key, value] of Object.entries(experiments)) {
                experimentsArray.push({ id: key, name: value });
            }
            state.proteinDesign.experiments = experimentsArray;
        },
        proteinDesign_inference(state, experiment) {
            state.proteinDesign.experiment = experiment;
        },
        proteinDesign_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.id = id;
            state.proteinDesign.experiments.push(experiment);
        },
        proteinDesign_deleteExperiment(state, experiment) {
            if (state.proteinDesign.experiment && state.proteinDesign.experiment.id === experiment.id) {
                state.proteinDesign.experiment = null;
            }
            state.proteinDesign.experiments = state.proteinDesign.experiments.filter(exp => exp.id !== experiment.id);
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
            formData.append('experimentId', experiment.id ?? '');
            formData.append('experimentName', experiment.name);
            const response = await axios({
                method: 'post',
                url: apiConstants.aminoAcid.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.aminoAcid.inference.mutation, response.data);
        },
        async aminoAcid_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.aminoAcid.experiments.path);
            commit(apiConstants.aminoAcid.experiments.mutation, response.data);
        },
        async aminoAcid_deleteExperiment({ commit }, { experiment }) {
            await axios.delete(apiConstants.aminoAcid.deleteExperiment.path, { data: { id: experiment.id } });
            commit(apiConstants.aminoAcid.deleteExperiment.mutation, experiment);
        },
        async aminoAcid_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperiment.path, { params: { id: experiment.id, name: experiment.name } });
            commit(apiConstants.aminoAcid.loadExperiment.mutation, { experiment: response.data });
        },
        async aminoAcid_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment)
                await axios.post(apiConstants.aminoAcid.changeExperimentName.path, { id: experiment.id, name: experiment.name });
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
        },
        async proteinDesign_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.proteinDesign.generateId.path);
            commit(apiConstants.proteinDesign.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async proteinDesign_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.id ?? '');
            formData.append('experimentName', experiment.name);
            const response = await axios({
                method: 'post',
                url: apiConstants.proteinDesign.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
            commit(apiConstants.proteinDesign.inference.mutation, response.data);
        },
        async proteinDesign_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.proteinDesign.experiments.path);
            commit(apiConstants.proteinDesign.experiments.mutation, response.data);
        },
        async proteinDesign_deleteExperiment({ commit }, { experiment }) {
            await axios.delete(apiConstants.proteinDesign.deleteExperiment.path, { data: { id: experiment.id } });
            commit(apiConstants.proteinDesign.deleteExperiment.mutation, experiment);
        },
        async proteinDesign_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.proteinDesign.loadExperiment.path, { params: { id: experiment.id, name: experiment.name } });
            commit(apiConstants.proteinDesign.loadExperiment.mutation, { experiment: response.data });
        },
        async proteinDesign_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment)
                await axios.post(apiConstants.proteinDesign.changeExperimentName.path, { id: experiment.id, name: experiment.name });
        },
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
    },
    proteinDesign: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.proteinDesign.experiments.action);
        },
        inference: async (payload) => {
            return await store.dispatch(apiConstants.proteinDesign.inference.action, { payload });
        },
        deleteExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.proteinDesign.deleteExperiment.action, { experiment });
        },
        loadExperiment: async (experiment) => {
            return await store.dispatch(apiConstants.proteinDesign.loadExperiment.action, { experiment });
        },
        addExperiment: () => {
            store.dispatch(apiConstants.proteinDesign.addExperiment.action, { experiment: { name: 'New experiment' } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.proteinDesign.changeExperimentName.action, { experiment });
        }
    }
}

export default store;