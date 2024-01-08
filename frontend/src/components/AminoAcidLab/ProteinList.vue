<script>
export default {
    props: ['experiment', 'api'],
    data() {
        return {
            proteinProgress: {},
            lastFetchedProgress: {}, // Store the last fetched progress data
            pollingInterval: null,
            isLoadingExperimentData: false // New property to track loading state
        };
    },
    methods: {
        loadPredictionsData(proteinId) {
            this.isLoadingExperimentData = true; // Set loading to true when starting to load data
            this.api.loadResults(this.experiment, proteinId).finally(() => {
                this.isLoadingExperimentData = false; // Set loading to false when data is loaded
            });
        },
        fetchAndUpdateProteinIds() {
            // Fetch the latest list of protein IDs from the API
            this.api.loadExperiment(this.experiment);
        },
        isProcessing(proteinId) {
            const progress = this.experiment.metaData.proteinIds[proteinId].progress.progress || 0;
            return progress < 100;
        },
        getButtonStyle(proteinId) {
            const progress = this.experiment.metaData.proteinIds[proteinId].progress.progress || 0;
            const baseStyle = {
                border: '2px solid green', 
                color: this.isNightMode ? 'white' : 'black',
                cursor: progress < 100 ? 'not-allowed' : 'pointer',
                margin: '5px'
            };
            
            if (progress < 100) {
                // Pulsing effect for ongoing progress
                return {
                    ...baseStyle,
                    border: '2px solid gray', 
                    color: this.isNightMode ? 'gray' : 'gray',
                };
            } else {
                return {
                    ...baseStyle,
                };
            }
        },
        startPollingProgress() {
            this.fetchAndUpdateProteinIds(); // Fetch immediately when component is mounted
            this.pollingInterval = setInterval(() => {
                this.fetchAndUpdateProteinIds(); // Continue fetching at regular intervals
            }, 5000);
        },
    },
    computed: {
        proteinIds() {
            return this.experiment.metaData.proteinIds ? Object.keys(this.experiment.metaData.proteinIds) : [];
        },
        isNightMode() {
            // Use window.matchMedia to check if the user prefers a dark color scheme
            return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
        },
        experimentEmpty() {
            if (!this.state.experiment.data)
                return true;

            if (typeof this.experiment.data === 'Array')
                return this.experiment.data.length === 0;

            return Object.keys(this.experiment.data).length === 0;
        },
    },
    mounted() {
        this.startPollingProgress();
    },
    beforeDestroy() {
        clearInterval(this.pollingInterval);
    }
}
</script>

<template>
    <div v-if="isLoadingExperimentData" class="loading-overlay">
            <div class="spinner-grow" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
            <p>Loading experiment data...</p>
    </div>

    <h4>Protein List:</h4>
    <div v-for="proteinId in proteinIds" :key="proteinId">
        <div class="container-fluid row align-items-center">
            <button class="btn"
                :style="getButtonStyle(proteinId)"
                :disabled="isLoadingExperimentData || isProcessing(proteinId)"
                @click="loadPredictionsData(proteinId)">
                <span v-if="isProcessing(proteinId)" class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>
                <span v-if="isProcessing(proteinId)">Processing... </span>
                {{ proteinId }}
            </button>
        </div>
    </div>  
</template>

<style>
.status-icon {
    margin-right: 10px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 24px;
}
</style>