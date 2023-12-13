<script>
export default {
    props: ['experiment', 'api'],
    data() {
        return {
            proteinProgress: {},
            lastFetchedProgress: {}, // Store the last fetched progress data
            pollingInterval: null
        };
    },
    methods: {
        loadPredictionsData(proteinId) {
            this.api.loadResults(this.experiment, proteinId);
        },
        fetchAndUpdateProteinIds() {
            // Fetch the latest list of protein IDs from the API
            debugger;
            this.api.loadExperiment(this.experiment.metaData);
        },
        isProcessing(proteinId) {
            const progress = this.experiment.metaData.proteinIds[proteinId].progress.progress || 0;
            return progress < 100;
        },
        getButtonStyle(proteinId) {
            const progress = this.experiment.metaData.proteinIds[proteinId].progress.progress || 0;
            if (progress < 100) {
                // Pulsing effect for ongoing progress
                return {
                    background: `linear-gradient(to right, green, transparent);`,
                    borderColor: 'green', // Optional: change the border color
                    color: 'white', // Optional: change the text color
                    cursor: 'not-allowed', // Change cursor to indicate it's not clickable
                    margin: '5px'
                };
            } else {
                // Static progress bar for completed progress
                return {
                    //background: `linear-gradient(to right, green ${progress}%, transparent ${progress}%)`,
                    background: `linear-gradient(to right, green, transparent);`,
                    borderColor: 'green', // Optional
                    color: 'white', // Optional
                    cursor: 'pointer', // Change cursor to indicate it's clickable
                    margin: '5px'
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
    <h4>Protein List:</h4>
    <div v-for="proteinId in proteinIds" :key="proteinId">
        <div class="container-fluid row align-items-center">
            <div class="status-icon">
                <div v-if="isProcessing(proteinId)" class="spinner"></div>
            </div>
            <button class="btn"
                :style="getButtonStyle(proteinId)"
                :disabled="isProcessing(proteinId)"
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

@keyframes spinner {
    to {transform: rotate(360deg);}
}

.spinner {
    width: 24px;
    height: 24px;
    border: 4px solid rgba(0, 0, 0, 0.1);
    border-radius: 50%;
    border-top-color: #007bff;
    animation: spinner 1s linear infinite;
}
</style>