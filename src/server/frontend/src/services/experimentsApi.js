import store from "../storage";

const ExperimentsApi = {
    getAllExperiments: async function () {
        const experimentsData1 = structuredClone(store.state.drugTargetData.inference);
        const experimentsData2 = structuredClone(store.state.drugTargetData.inference);
        experimentsData1.name = 'Test1';
        experimentsData2.name = 'Test2';
        return [experimentsData1, experimentsData2];
    },
    deleteExperiment: async function (experiment) {
        return experiments.splice(experiments.indexOf(experiment), 1);
    },
    saveExperiment: async function (experiment) {
        const newExperiment = structuredClone(drugTargetData.inference);
        newExperiment.name = 'New experiment';
        store.state.drugTargetData.experiments.push(newExperiment);
    }
}

export default ExperimentsApi;