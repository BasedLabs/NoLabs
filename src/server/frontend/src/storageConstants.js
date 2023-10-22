const baseUrl = 'http://127.0.0.1:5000'

export const apiConstants = {
    aminoAcid: {
        addExperiment: {
            mutation: 'aminoAcid_addExperiment',
            action: 'aminoAcid_addExperiment'
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
            mutation: 'aminoAcid_eleteExperiment',
            action: 'aminoAcid_deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + '/api/amino-acid/load-experiment',
            mutation: 'aminoAcid_loadExperiment',
            action: 'aminoAcid_loadExperiment'
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
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/drug-target/experiments',
            mutation: 'drugTarget_getAllExperiments',
            action: 'drugTarget_getAllExperiments'
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
        changeExperimentName: {
            path: baseUrl + '/api/drug-target/change-experiment-name',
            action: 'drugTarget_changeExperimentName'
        }
    }
}

export default apiConstants;