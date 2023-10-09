const baseUrl = 'http://127.0.0.1:5000'

export const apiConstants = {
    aminoAcid: {
        inference: {
            path: baseUrl + '/api/amino-acid-inference',
            mutation: 'aminoAcidInference',
            action: 'aminoAcidInference'
        }
    },
    drugTargetDiscovery: {
        experiments: {
            path: baseUrl + '/api/drug-target-discovery-experiments',
            mutation: 'getAllExperiments',
            action: 'getAllExperiments'
        },
        addExperiment: {
            path: baseUrl + 'api/drug-target-discovery-add-experiment',
            mutation: 'addExperiment',
            action: 'addExperiment'
        },
        deleteExperiment: {
            path: baseUrl + 'api/drug-target-discovery-delete-experiment',
            mutation: 'deleteExperiment',
            action: 'deleteExperiment'
        },
        loadExperiment: {
            path: baseUrl + 'api/drug-target-discovery-load-experiment',
            mutation: 'loadExperiment',
            action: 'loadExperiment'
        }
    }
}

export default apiConstants;