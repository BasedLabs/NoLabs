<script>
import LigandViewer from './LigandViewer.vue' // Assuming you have a similar viewer for ligands

export default {
    props: ['experiment', 'api'],
    components: {
        LigandViewer
    },
    data() {
        return {
            isLoadingExperimentData: false,
            ligands: {},
            isModalOpen: false,
            selectedLigand: null
        };
    },
    methods: {
        loadLigands() {
            this.api.loadLigands(this.experiment)
                .then(() => {
                    this.ligands = this.experiment.ligands;
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
                    if (file.name.endsWith('.sdf')) {
                        this.uploadLigandFile(file);
                    } else {
                        console.error('Invalid file type:', file.name);
                    }
                });
            }
        },
        uploadLigandFile(file) {
            this.api.addLigand(this.experiment, file)
                .then(() => {
                    this.loadLigands();
                })
                .catch(error => {
                    console.error('Error uploading file:', error);
                });
        },
        deleteLigand(ligandId) {
            this.api.deleteLigand(this.experiment.metaData.id, ligandId)
                .then(() => {
                    // Refresh the targets list after deletion
                    this.loadLigands();
                })
                .catch(error => {
                    console.error('Error deleting ligand:', error);
                });
        },
        openModal(ligand) {
            this.selectedLigand = ligand;
            this.isModalOpen = true;
        },
        closeModal() {
            this.isModalOpen = false;
        }
    },
    mounted() {
        this.loadLigands();
    }
}
</script>

<template>
    <div>
        <!-- Loading Indicator -->
        <div v-if="isLoadingExperimentData" class="loading-overlay">
            <div class="spinner-grow" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
            <p>Loading experiment data...</p>
        </div>
        
        <div class="file-upload-area" @dragover.prevent="dragOverHandler" @drop.prevent="dropHandler">
            <p>Drag and drop your .sdf files here, or click to select files</p>
            <input type="file" id="fileInput" multiple @change="handleFileUpload" accept=".sdf" style="display: none;" />
            <label for="fileInput" class="file-upload-button">Select Files</label>
        </div>

        <!-- Display Uploaded Ligands -->
        <div class="container-fluid" v-if="Object.keys(ligands).length > 0">
            <h4>Uploaded ligands: </h4>
            <div v-for="(ligand, id) in ligands" :key="id" class="ligand-button">
                <button class="btn btn-info" style="width: 50vw;" @click="openModal(ligand)">{{ ligand.metadata.name }}</button>
                <button class="btn btn-danger" @click="deleteLigand(id)">Delete</button>
            </div>
        </div>

        <!-- Ligand Viewer Modal -->
        <div v-if="isModalOpen" class="modal">
            <LigandViewer :api="this.api" :experiment="this.experiment" :ligand="selectedLigand" @close="closeModal" />
        </div>
    </div>
</template>
