export const baseUrl = 'http://localhost:5000'

export const apiConstants = {
    conformations: {
        addExperiment: {
            mutation: 'conformations_addExperiment',
            action: 'conformations_addExperiment'
        },
        generateId: {
            path: baseUrl + '/api/conformations/generate-id'
        },
        experiments: {
            path: baseUrl + '/api/conformations/experiments',
            mutation: 'conformations_getAllExperiments',
            action: 'conformations_getAllExperiments'
        },
        inference: {
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
}

export default apiConstants;