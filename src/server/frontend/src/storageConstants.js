const baseUrl = 'http://127.0.0.1:5000'

export const apiConstants = {
    aminoAcid: {
        inference: {
            path: baseUrl + '/api/amino-acid/inference',
            mutation: 'aminoAcidInference',
            action: 'aminoAcidInference'
        }
    },
    // When we click on submit - we trigger server inference, then we save it, and then we return this inference to the UI
    // inference model is a complete experiment model
    // getExperiment - get the experiment
    // experiments - short data with id and name fields
    // delete experiment - pass id and delete
    // experiment name - is Experiment # - generated

    drugTargetDiscovery: {
        experiments: { // response: [{name: 'exp1'}, {name: 'exp2'}]
            path: baseUrl + '/api/drug-target/experiments',
            mutation: 'getAllExperiments',
            action: 'getAllExperiments'
        },
        inference: { // request: {name: 'exp1', sdf: 'sdf', pdb: 'pdb'} response: {exp, data: []}
            path: baseUrl + 'api/drug-target/inference',
            mutation: 'drugTargetInference',
            action: 'drugTargetInference'
        },
        deleteExperiment: { // request: {name: '123'}
            path: baseUrl + 'api/drug-target/delete-experiment',
            mutation: 'deleteExperiment',
            action: 'deleteExperiment'
        },
        loadExperiment: { // request: {name: '123'}, response: {exp, data: []}
            path: baseUrl + 'api/drug-target/load-experiment',
            mutation: 'loadExperiment',
            action: 'loadExperiment'
        }
    }
}

export default apiConstants;