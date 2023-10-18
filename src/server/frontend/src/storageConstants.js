const baseUrl = 'http://127.0.0.1:5000'

export const apiConstants = {
    aminoAcid: {
        addExperiment: {
            mutation: 'aminoAcid.addExperiment',
            action: 'aminoAcid.addExperiment'
        },
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/amino-acid/experiments',
            mutation: 'aminoAcid.getAllExperiments',
            action: 'aminoAcid.getAllExperiments'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + '/api/amino-acid/inference',
            mutation: 'aminoAcid.inference',
            action: 'aminoAcid.inference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + '/api/amino-acid/delete-experiment',
            mutation: 'daminoAcid.eleteExperiment',
            action: 'aminoAcid.deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + '/api/amino-acid/load-experiment',
            mutation: 'aminoAcid.loadExperiment',
            action: 'aminoAcid.loadExperiment'
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
            mutation: 'drugTarget.addExperiment',
            action: 'drugTarget.addExperiment'
        },
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/drug-target/experiments',
            mutation: 'drugTarget.getAllExperiments',
            action: 'drugTarget.getAllExperiments'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + '/api/drug-target/inference',
            mutation: 'drugTarget.inference',
            action: 'drugTarget.inference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + '/api/drug-target/delete-experiment',
            mutation: 'drugTarget.deleteExperiment',
            action: 'drugTarget.deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + '/api/drug-target/load-experiment',
            mutation: 'drugTarget.loadExperiment',
            action: 'drugTarget.loadExperiment'
        }
    }
}

export default apiConstants;