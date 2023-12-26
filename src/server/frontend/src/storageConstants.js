export const baseUrl = 'http://127.0.0.1:5000'

export const apiConstants = {
    aminoAcid: {
        addExperiment: {
            mutation: 'aminoAcid_addExperiment',
            action: 'aminoAcid_addExperiment'
        },
        generateId: {
            path: baseUrl + '/api/amino-acid/generate-id'
        },
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/amino-acid/experiments',
            mutation: 'aminoAcid_getAllExperiments',
            action: 'aminoAcid_getAllExperiments'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + '/api/amino-acid/inference',
            mutation: 'aminoAcid_inference',
            action: 'aminoAcid_inference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + '/api/amino-acid/delete-experiment',
            mutation: 'aminoAcid_deleteExperiment',
            action: 'aminoAcid_deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, protein_ids: []}
            path: baseUrl + '/api/amino-acid/load-experiment',
            mutation: 'aminoAcid_loadExperiment',
            action: 'aminoAcid_loadExperiment'
        },
        loadResults : {
            path: baseUrl + '/api/amino-acid/load-results',
            mutation: 'aminoAcid_loadResults',
            action: 'aminoAcid_loadResults'
        },
        loadExperimentProgress : {
            path: baseUrl + '/api/amino-acid/load-experiment-progress',
            mutation: 'aminoAcid_loadExperimentProgress',
            action: 'aminoAcid_loadExperimentProgress'
        },
        loadExperimentInstanceProgress : {
            path: baseUrl + '/api/amino-acid/load-experiment-instance-progress',
            mutation: 'aminoAcid_loadExperimentInstanceProgress',
            action: 'aminoAcid_loadExperimentInstanceProgress'
        },
        changeExperimentName: {
            path: baseUrl + '/api/amino-acid/change-experiment-name',
            action: 'aminoAcid_changeExperimentName'
        }
    },
    // When we click on submit - we trigger server inference, then we save it, and then we return this inference to the UI
    // inference model is a complete experiment model
    // getExperiment - get the experiment
    // experiments - short data with id and name fields
    // delete experiment - pass id and delete
    // experiment name - is Experiment # - generated

    drugTarget: {
        addExperiment: {
            mutation: 'drugTarget_addExperiment',
            action: 'drugTarget_addExperiment'
        },
        generateId: {
            path: baseUrl + '/api/drug-target/generate-id'
        },
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/drug-target/experiments',
            mutation: 'drugTarget_getAllExperiments',
            action: 'drugTarget_getAllExperiments'
        },
        addTarget: {
            path: baseUrl + '/api/drug-target/add-target',
            mutation: 'drugTarget_addTarget',
            action: 'drugTarget_addTarget'
        },
        loadTargets: {
            path: baseUrl + '/api/drug-target/load-targets',
            mutation: 'drugTarget_loadTargets',
            action: 'drugTarget_loadTargets'
        },
        setBindingPocket: {
            path: baseUrl + '/api/drug-target/set-binding-pocket',
            mutation: 'drugTarget_setBindingPocket',
            action: 'drugTarget_setBindingPocket'
        },
        loadBindingPocket: {
            path: baseUrl + '/api/drug-target/get-binding-pocket'
        },
        predictBindingPocket: {
            path: baseUrl + '/api/drug-target/predict-binding-pocket'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + '/api/drug-target/inference',
            mutation: 'drugTarget_inference',
            action: 'drugTarget_inference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + '/api/drug-target/delete-experiment',
            mutation: 'drugTarget_deleteExperiment',
            action: 'drugTarget_deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + '/api/drug-target/load-experiment',
            mutation: 'drugTarget_loadExperiment',
            action: 'drugTarget_loadExperiment'
        },
        loadResults : {
            path: baseUrl + '/api/drug-target/load-results',
            mutation: 'drugTarget_loadResults',
            action: 'drugTarget_loadResults'
        },
        changeExperimentName: {
            path: baseUrl + '/api/drug-target/change-experiment-name',
            action: 'drugTarget_changeExperimentName'
        },
        downloadCombinedPdb: {
            path: baseUrl + '/api/drug-target/download-combined-pdb'
        }
    },

    conformations: {
        addExperiment: {
            mutation: 'conformations_addExperiment',
            action: 'conformations_addExperiment'
        },
        generateId: {
            path: baseUrl + '/api/conformations/generate-id'
        },
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/conformations/experiments',
            mutation: 'conformations_getAllExperiments',
            action: 'conformations_getAllExperiments'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + '/api/conformations/inference',
            mutation: 'conformations_inference',
            action: 'conformations_inference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + '/api/conformations/delete-experiment',
            mutation: 'conformations_deleteExperiment',
            action: 'conformations_deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + '/api/conformations/load-experiment',
            mutation: 'conformations_loadExperiment',
            action: 'conformations_loadExperiment'
        },
        changeExperimentName: {
            path: baseUrl + '/api/conformations/change-experiment-name',
            action: 'conformations_changeExperimentName'
        }
    },

    proteinDesign: {
        addExperiment: {
            mutation: 'proteinDesign_addExperiment',
            action: 'proteinDesign_addExperiment'
        },
        generateId: {
            path: baseUrl + '/api/protein-design/generate-id'
        },
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/protein-design/experiments',
            mutation: 'proteinDesign_getAllExperiments',
            action: 'proteinDesign_getAllExperiments'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + '/api/protein-design/inference',
            mutation: 'proteinDesign_inference',
            action: 'proteinDesign_inference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + '/api/protein-design/delete-experiment',
            mutation: 'proteinDesign_deleteExperiment',
            action: 'proteinDesign_deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + '/api/protein-design/load-experiment',
            mutation: 'proteinDesign_loadExperiment',
            action: 'proteinDesign_loadExperiment'
        },
        changeExperimentName: {
            path: baseUrl + '/api/protein-design/change-experiment-name',
            action: 'proteinDesign_changeExperimentName'
        }
    }
}

export default apiConstants;