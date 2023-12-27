import { createStore } from 'vuex';
import axios from 'axios';
import apiConstants from './storageConstants';

const store = createStore({
    state: {
        aminoAcid: {
            experiment: {
                metaData: null,
                data: null
            },
            experiments: []
        },
        drugTarget: {
            experiment: {
                metaData: null,
                targets: null,
                data: null
            },
            experiments: []
        },
        conformations: {
            experiment: {
                metaData: null,
                data: null
            },
            experiments: []
        },
        proteinDesign: {
            experiment: {
                metaData: null,
                data: null,
                errors: []
            },
            experiments: []
        }
    },
    mutations: {
        aminoAcid_loadExperiment(state, { experiment }) {
            state.aminoAcid.experiment.metaData = experiment;
        },
        aminoAcid_loadResults(state, { experimentData }) {
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
                    metaData: {
                        id: key,
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
            state.aminoAcid.experiments = state.aminoAcid.experiments.filter(exp => exp.metaData.id !== experiment.metaData.id);
        },
        drugTarget_loadExperiment(state, { experiment }) {
            //state.drugTarget.experiment = experiment;
            state.drugTarget.experiment.metaData = experiment;
        },
        drugTarget_loadResults(state, { experimentData }) {
            state.drugTarget.experiment.data = experimentData.data;
            state.drugTarget.experiment.metaData['selectedLigandId'] = experimentData.ligandId;
            state.drugTarget.experiment.metaData['selectedProteinId'] = experimentData.proteinId;
        },
        drugTarget_getAllExperiments(state, experiments) {
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
            state.drugTarget.experiments = experimentsArray;
        },
        drugTarget_loadTargets(state, { data }) {
            //state.drugTarget.experiment = experiment;
            state.drugTarget.experiment.targets = data.targets;
        },
        drugTarget_inference(state, experiment) {
            state.drugTarget.experiment = experiment;
        },
        drugTarget_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.metaData.id = id;
            state.drugTarget.experiments.push(experiment);
        },
        drugTarget_deleteExperiment(state, experiment) {
            if (state.drugTarget.experiment && state.drugTarget.experiment.metaData.id === experiment.metaData.id) {
                state.drugTarget.experiment = null;
            }
            state.drugTarget.experiments = state.drugTarget.experiments.filter(exp => exp.metaData.id !== experiment.metaData.id);
        },
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
            state.conformations.experiments = state.conformations.experiments.filter(exp => exp.id !== experiment.metaData.id);
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
                experimentsArray.push({
                    metaData: {
                        id: key,
                        name: value.name,
                        progress: value.progress,
                        date: value.date
                    }
                });
            }
            state.proteinDesign.experiments = experimentsArray;
        },
        proteinDesign_inference(state, experiment) {
            state.proteinDesign.experiment = experiment;
        },
        proteinDesign_addExperiment(state, data) {
            const { experiment, id } = data;
            experiment.metaData.id = id;
            state.proteinDesign.experiments.push(experiment);
        },
        proteinDesign_deleteExperiment(state, experiment) {
            state.proteinDesign.experiments = state.proteinDesign.experiments.filter(exp => exp.metaData.id !== experiment.metaData.id);
        }
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
                {
                    params:
                    {
                        id: experiment.metaData.id,
                        name: experiment.metaData.name,
                        proteinId: proteinId
                    }
                });
            commit(apiConstants.aminoAcid.loadResults.mutation, { experimentData: response.data });
        },
        async aminoAcid_loadExperimentProgress({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperimentProgress.path,
                { params: { id: experiment.metaData.id } });
            commit(apiConstants.aminoAcid.loadExperimentProgress.mutation, { experimentData: response.data });
        },
        async aminoAcid_loadExperimentInstanceProgress({ commit }, { experiment, proteinId }) {
            const response = await axios.get(apiConstants.aminoAcid.loadExperimentInstanceProgress.path,
                {
                    params:
                    {
                        id: experiment.metaData.id,
                        proteinId: proteinId
                    }
                });
            commit(apiConstants.aminoAcid.loadExperimentInstanceProgress.mutation, { experimentInstanceProgress: response.data });
        },
        async aminoAcid_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment.metaData)
                await axios.post(apiConstants.aminoAcid.changeExperimentName.path, { id: experiment.metaData.id, name: experiment.metaData.name });
        },
        async drugTarget_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.drugTarget.generateId.path);
            commit(apiConstants.drugTarget.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async drugTarget_addTarget({ commit }, { experiment, file }) {
            // Constructing form data for file upload
            let formData = new FormData();
            formData.append('proteinFileInput', file);
            formData.append('experimentId', experiment.metaData.id);
            formData.append('experimentName', experiment.metaData.name);
            // Assuming there's an endpoint in your API for uploading target files
            return await axios({
                method: 'post',
                url: apiConstants.drugTarget.addTarget.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
        },
        async drugTarget_loadTargets({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.drugTarget.loadTargets.path, { params: { id: experiment.metaData.id, name: experiment.metaData.name } });
            commit(apiConstants.drugTarget.loadTargets.mutation, { data: response.data });
        },
        async drugTarget_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.metaData.id ?? '');
            formData.append('experimentName', experiment.metaData.name);
            const response = await axios({
                method: 'post',
                url: apiConstants.drugTarget.inference.path,
                data: formData,
                headers: { "Content-Type": "multipart/form-data" },
            });
        },
        async drugTarget_getAllExperiments({ commit }) {
            const response = await axios.get(apiConstants.drugTarget.experiments.path);
            commit(apiConstants.drugTarget.experiments.mutation, response.data);
        },
        async drugTarget_deleteExperiment({ commit }, { experiment }) {
            await axios.delete(apiConstants.drugTarget.deleteExperiment.path, { data: { id: experiment.metaData.id } });
            commit(apiConstants.drugTarget.deleteExperiment.mutation, experiment);
        },
        async drugTarget_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.drugTarget.loadExperiment.path, { params: { id: experiment.metaData.id, name: experiment.metaData.name } });
            commit(apiConstants.drugTarget.loadExperiment.mutation, { experiment: response.data });
        },
        async drugTarget_loadResults({ commit }, { experiment, proteinId, ligandId }) {
            const response = await axios.get(apiConstants.drugTarget.loadResults.path,
                {
                    params:
                    {
                        id: experiment.metaData.id,
                        name: experiment.metaData.name,
                        proteinId: proteinId,
                        ligandId: ligandId
                    }
                });
            commit(apiConstants.drugTarget.loadResults.mutation, { experimentData: response.data });
        },
        async drugTarget_changeExperimentName({ commit }, { experiment }) {
            if ('id' in experiment.metaData)
                await axios.post(apiConstants.drugTarget.changeExperimentName.path, { id: experiment.metaData.id, name: experiment.metaData.name });
        },
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
        },
        async proteinDesign_addExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.proteinDesign.generateId.path);
            commit(apiConstants.proteinDesign.addExperiment.mutation, { experiment, id: response.data.id });
        },
        async proteinDesign_inference({ commit }, { payload }) {
            const { form, experiment } = payload;
            const formData = new FormData(form);
            formData.append('experimentId', experiment.metaData.id);
            formData.append('experimentName', experiment.metaData.name);
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
            await axios.delete(apiConstants.proteinDesign.deleteExperiment.path, { data: { id: experiment.metaData.id } });
            commit(apiConstants.proteinDesign.deleteExperiment.mutation, experiment);
        },
        async proteinDesign_loadExperiment({ commit }, { experiment }) {
            const response = await axios.get(apiConstants.proteinDesign.loadExperiment.path, { params: { id: experiment.metaData.id, name: experiment.metaData.name } });
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
            store.dispatch(apiConstants.aminoAcid.addExperiment.action, { experiment: { metaData: { name: `New experiment` } } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.aminoAcid.changeExperimentName.action, { experiment });
        }
    },
    drugTarget: {
        getAllExperiments: async () => {
            return await store.dispatch(apiConstants.drugTarget.experiments.action);
        },
        addTarget: async (experiment, file) => {
            return await store.dispatch(apiConstants.drugTarget.addTarget.action, { experiment, file });
        },
        loadTargets: async (experiment) => {
            return await store.dispatch(apiConstants.drugTarget.loadTargets.action, { experiment });
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
        loadResults: async (experiment, proteinId, ligandId) => {
            return await store.dispatch(apiConstants.drugTarget.loadResults.action, { experiment, proteinId, ligandId });
        },
        addExperiment: () => {
            store.dispatch(apiConstants.drugTarget.addExperiment.action, { experiment: { metaData: { name: `New experiment` } } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.drugTarget.changeExperimentName.action, { experiment });
        },
        downloadCombinedPdb: async (experiment, selectedProteinId, selectedLigandId) => {
            return await axios.post(apiConstants.drugTarget.downloadCombinedPdb.path,
                { experimentId: experiment.metaData.id, proteinId: selectedProteinId, ligandId: selectedLigandId });
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
            store.dispatch(apiConstants.conformations.addExperiment.action, { experiment: { metaData: { name: `New experiment` } } });
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
            store.dispatch(apiConstants.proteinDesign.addExperiment.action, { experiment: { metaData: { name: `New experiment` } } });
        },
        changeExperimentName: async (experiment) => {
            await store.dispatch(apiConstants.proteinDesign.changeExperimentName.action, { experiment });
        }
    }
}

export default store;