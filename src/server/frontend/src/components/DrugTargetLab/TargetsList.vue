<script>
import TargetProteinViewer from './TargetProteinViewer.vue'

export default {
    props: ['experiment', 'api'],
    components: {
        TargetProteinViewer
    },
    data() {
        return {
            proteinProgress: {},
            lastFetchedProgress: {},
            pollingInterval: null,
            isLoadingExperimentData: false,
            targets: {},
            isModalOpen: false,
            selectedTarget: null
        };
    },
    methods: {
        loadTargets() {
            this.api.loadTargets(this.experiment)
                .then(() => {
                    this.targets = this.experiment.targets;// Update the targets list with the new file
                    // Optionally, refresh targets from backend if needed
                    // this.loadTargets(); 
                });
        },
        dragOverHandler(event) {
            event.dataTransfer.dropEffect = 'copy';
        },
        dropHandler(event) {
            this.handleFileUpload(event);
        },
        handleFileUpload(event) {
            const files = event.dataTransfer ? event.dataTransfer.files : event.target.files;
            if (files && files.length > 0) {
                Array.from(files).forEach(file => {
                    if (file.name.endsWith('.fasta')) {
                        this.uploadProteinFile(file);
                    } else {
                        console.error('Invalid file type:', file.name);
                    }
                });
            }
        },
        uploadProteinFile(file) {
            this.api.addTarget(this.experiment, file)
                .then(() => {
                    this.loadTargets();
                })
                .catch(error => {
                    console.error('Error uploading file:', error);
                });
        },
        deleteTarget(targetId) {
            this.api.deleteTarget(this.experiment.metaData.id, targetId)
                .then(() => {
                    // Refresh the targets list after deletion
                    this.loadTargets();
                })
                .catch(error => {
                    console.error('Error deleting target:', error);
                });
        },
        openModal(target) {
            this.selectedTarget = target;
            this.isModalOpen = true;
        },
        closeModal() {
            this.isModalOpen = false;
        },
    },
    computed: {
        isNightMode() {
            return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
        }
    },
    mounted() {
        this.loadTargets();
    },
    beforeDestroy() {
        clearInterval(this.pollingInterval);
    }
    }
</script>

<template>
    <div>
        <div v-if="isLoadingExperimentData" class="loading-overlay">
            <div class="spinner-grow" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
            <p>Loading experiment data...</p>
        </div>

        <div class="file-upload-area" @dragover.prevent="dragOverHandler" @drop.prevent="dropHandler">
            <p>Drag and drop your .fasta files here, or click to select files</p>
            <input type="file" id="fileInput" multiple @change="handleFileUpload" accept=".fasta" style="display: none;" />
            <label for="fileInput" class="file-upload-button">Select Files</label>
        </div>

        <div class="container-fluid" v-if="Object.keys(targets).length > 0">
            <h4>Uploaded targets: </h4>
            <div v-for="(target, id) in targets" :key="id" class="target-button">
                <button class="btn btn-info" style="width: 50vw;" @click="openModal(target)">{{ target.metadata.name }}</button>
                <button class="btn btn-danger" @click="deleteTarget(id)">Delete</button>
            </div>
        </div>

        <!-- Modal -->
        <div v-if="isModalOpen" class="modal">
            <TargetProteinViewer :api="this.api" :experiment="this.experiment" :target="selectedTarget" @close="closeModal" />
        </div>
    </div>
</template>

<style>

.file-upload-area {
    border: 2px dashed #007bff;
    border-radius: 5px;
    text-align: center;
    padding: 20px;
    margin: 10px 0;
    cursor: pointer;
    background-color: #3e3434;
    transition: background-color 0.3s ease-in-out;
}

.file-upload-area:hover {
    background-color: #494b4d;
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