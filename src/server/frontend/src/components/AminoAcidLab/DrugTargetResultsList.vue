<script>
import ResultProteinViewer from './ResultProteinViewer.vue'


export default {
    props: ['experiment', 'api'],
    components: {
        ResultProteinViewer
    },
    data() {
        return {
            proteinProgress: {},
            lastFetchedProgress: {}, // Store the last fetched progress data
            dropdownOpen: {}, // Object to track open state of each dropdown
            pollingInterval: null,
            isLoadingExperimentData: false // New property to track loading state
        };
    },
    methods: {
        loadPredictionsData(proteinId, ligandId) {
            this.isLoadingExperimentData = true; // Set loading to true when starting to load data
            this.api.loadPredictionData(this.experiment, proteinId, ligandId).finally(() => {
                this.isLoadingExperimentData = false; // Set loading to false when data is loaded
                this.isModalOpen = true;
            });
        },
        fetchAndUpdateProteinIds() {
            // Fetch the latest list of protein IDs from the API
            this.api.loadResults(this.experiment);
        },
        toggleDropdown(proteinId) {
            if (this.dropdownOpen[proteinId] === undefined) {
            // Initialize as open if not already present
            this.dropdownOpen[proteinId] = true;
            } else {
                // Toggle the current state
                this.dropdownOpen[proteinId] = !this.dropdownOpen[proteinId];
            }
        },
        getButtonStyle() {
            const baseStyle = {
                border: '2px solid green', 
                color: this.isNightMode ? 'white' : 'black',
                cursor: 'not-allowed',
                margin: '5px'
            };
            return {
                ...baseStyle,
                border: '2px solid gray', 
                color: this.isNightMode ? 'gray' : 'gray',
            };
        },
        startPollingProgress() {
            this.fetchAndUpdateProteinIds(); // Fetch immediately when component is mounted
            this.pollingInterval = setInterval(() => {
                this.fetchAndUpdateProteinIds(); // Continue fetching at regular intervals
            }, 5000);
        },
        ligandIds(proteinId) {
            if (this.experiment.results && this.experiment.results.proteinIds && this.experiment.results.proteinIds[proteinId]) {
                return this.experiment.results.proteinIds[proteinId].ligandIds 
                    ? Object.values(this.experiment.results.proteinIds[proteinId].ligandIds) 
                    : [];
            }
            return [];
        },
        closeModal() {
            this.isModalOpen = false;
        },
        isTargetInQueue(proteinId) {
            const targetLigands = this.experiment.ligands 
                ? Object.keys(this.experiment.ligands)
                : [];
            const resultLigands = this.experiment.results && this.experiment.results.proteinIds && this.experiment.results.proteinIds[proteinId] && this.experiment.results.proteinIds[proteinId].ligandIds
                ? Object.keys(this.experiment.results.proteinIds[proteinId].ligandIds)
                : [];

            return !this.arraysEqual(targetLigands.sort(), resultLigands.sort());
        },
        arraysEqual(a, b) {
            if (a.length !== b.length) return false;
            return a.every((val) => b.includes(val));
        },
        isInResultsOrProcessing(targetId, ligandId) {
            // Check if target-ligand pair is in results or currently processing
            return this.experiment.results 
                && this.experiment.results.proteinIds 
                && this.experiment.results.proteinIds[targetId]
                && (this.experiment.results.proteinIds[targetId].ligandIds.includes(ligandId)
                    || (this.experiment.results.proteinIds[targetId].ligandResultsAvailable
                        && this.experiment.results.proteinIds[targetId].ligandResultsAvailable[ligandId] === false));
        },
        isProcessing(proteinId) {
            // Check if any ligand for the given proteinId is still processing
            return this.ligandIds(proteinId).some(ligandId => !this.isResultAvailable(proteinId, ligandId));
        },
        isResultAvailable(proteinId, ligandId) {
            if (this.experiment.results.proteinIds
            && this.experiment.results.proteinIds[proteinId]
            && this.experiment.results.proteinIds[proteinId].ligandIds) {
            
            const ligandIndex = this.experiment.results.proteinIds[proteinId].ligandIds.indexOf(ligandId);
            
            if (ligandIndex !== -1 && this.experiment.results.proteinIds[proteinId].ligandResultsAvailable) {
                    return this.experiment.results.proteinIds[proteinId].ligandResultsAvailable[ligandIndex];
                }
            }
            return false; // Return false if the ligandId is not found or other conditions are not met
        },
        
        getLigandsForTarget(targetId) {
            // Return the ligands associated with the given target ID
            return this.experiment.ligands || [];
        },
    },
    computed: {
        targetsLoading() {
            const targetKeys = this.experiment && this.experiment.targets ? Object.keys(this.experiment.targets) : [];
            const resultKeys = this.proteinIds;
            return targetKeys.length !== resultKeys.length;
        },
        targetList() {
            return this.experiment && this.experiment.targets ? Object.keys(this.experiment.targets) : [];
        },
        proteinIds() {
            return this.experiment && this.experiment.results && this.experiment.results.proteinIds 
               ? Object.keys(this.experiment.results.proteinIds) 
               : [];
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

    <div>
        <div v-for="targetId in targetList" :key="targetId" class="container-fluid row align-items-center">
            <div v-for="ligandId in getLigandsForTarget(targetId)" :key="ligandId">
                <button v-if="!isInResultsOrProcessing(targetId, ligandId.metadata.name)" class="btn" :style="getButtonStyle()" :disabled="true">
                    <span>Queue > </span>
                    {{ this.experiment.targets[targetId].metadata.name }}: {{ ligandId.metadata.name }}
                </button>
            </div>
        </div>
    </div>

    <div v-for="proteinId in proteinIds" :key="proteinId">
        <div  v-if="isProcessing(proteinId)" class="container-fluid row align-items-center">
            <div v-for="ligandId in ligandIds(proteinId)" :key="ligandId">
                <button v-if="!isResultAvailable(proteinId, ligandId)" class="btn" :style="getButtonStyle()" :disabled="true">
                    <span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>
                    {{ this.experiment.targets[proteinId].metadata.name }}: {{ ligandId }}
                </button>
            </div>
        </div>
    </div>

    <h4 style="margin-top: 20px;">Results:</h4>
    <div v-for="proteinId in proteinIds" :key="proteinId">
        <div class="container-fluid row align-items-center">
            <div class="btn-group">
                <button class="btn btn-primary dropdown-toggle"
                        type="button" 
                        @click="toggleDropdown(proteinId)"
                        aria-expanded=false>
                    {{ this.experiment.targets[proteinId].metadata.name }}
                </button>
            </div>
            <div v-if="dropdownOpen[proteinId]" v-for="ligandId in ligandIds(proteinId)" :key="ligandId">
                    <button class="btn btn-primary"
                            @click="loadPredictionsData(proteinId, ligandId)">
                        {{ ligandId }}
                    </button>
            </div>
        </div>
    </div>

    <div v-if="isModalOpen" class="modal">
        <ResultProteinViewer :api="this.api" :experiment="this.experiment" @close="closeModal" />
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

.modal {
  position: fixed;
  left: 0;
  top: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(0, 0, 0, 0.5);
  display: flex;
  justify-content: center;
  align-items: center;
}
</style>