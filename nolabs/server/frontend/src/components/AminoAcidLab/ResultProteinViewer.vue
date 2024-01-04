<script>
import axios from 'axios';
export default {
    props: ['api', 'experiment'],
    data() {
        const views = {
            default: { key: 'default', title: 'Default representation' },
            cartoon: { key: 'cartoon', title: 'Cartoon' },
            backbone: { key: 'backbone', title: 'Backbone' },
            ballsAndSticks: { key: 'ball+stick', title: 'Balls and sticks' },
            contact: { key: 'contact', title: 'Contact' },
            helixorient: { key: 'helixorient', title: 'Helixorient' },
            hyperball: { key: 'hyperball', title: 'Hyperball' },
            licorice: { key: 'licorice', title: 'Licorice' },
            ribbon: { key: 'ribbon', title: 'Ribbon' },
            rope: { key: 'rope', title: 'Rope' },
            surface: { key: 'surface', title: 'Surface' },
            spacefill: { key: 'spacefill', title: 'Spacefill' },
            unitcell: { key: 'unitcell', title: 'Unitcell' }
        };
        return {
            views: views,
            viewsItems: Object.keys(views),
            selectedView: views.default,
            stage: null,
            pdbComponent: null,
            ligandComponent: null,
            selectedResidues: new Set(),
            sequence: [],
            selectionMode: false,
            showNoPocketModal: false,
            isLoading: false,
        }
    },
    methods: {
        setView(evt) {
            const viewKey = evt.target.value;
            this.selectedView = this.views[viewKey];
            if (this.selectedView === this.views.default) {
                this.render();
                return;
            }
            this.pdbComponent.removeAllRepresentations();
            this.pdbComponent.addRepresentation(this.selectedView.key);
        },
        cleanStage() {
            if (this.stage)
                this.stage.removeAllComponents();
        },
        render() {
            setTimeout(() => {
                // Render 3d structure
                document.getElementById('viewport').innerHTML = '';
                this.cleanStage();
                this.stage = new NGL.Stage("viewport");
                this.stage.setParameters({ backgroundColor: 'white',  tooltip: true });
                const proteinFileContentBlob = new Blob([this.experiment.data.pdb], { type: 'text/plain' });
                const ligandFileContentBlob = new Blob([this.experiment.data.sdf], { type: 'text/plain' });
                const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
                const sdfFile = new File([ligandFileContentBlob], 'ligand.sdf', { type: 'text/plain' });
                this.stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
                    this.pdbComponent = component;
                    this.stage.loadFile(sdfFile, { defaultRepresentation: true }).then((sdfComponent) => {
                        this.ligandComponent = sdfComponent;
                    });
                });
            }, 200);
        },
        downloadPdbFile(evt) {
            const filename = this.experiment.data.proteinName + "_" + this.experiment.data.ligandName  + '.pdb';
            const blob = new Blob([this.experiment.data.pdb], { type: 'text/plain' });
            if (window.navigator.msSaveOrOpenBlob) {
                window.navigator.msSaveBlob(blob, filename);
            } else {
                const elem = window.document.createElement('a');
                elem.href = window.URL.createObjectURL(blob);
                elem.download = filename;
                document.body.appendChild(elem);
                elem.click();
                document.body.removeChild(elem);
            }
        }
    },
    computed: {
    },
    mounted() {
        this.render();
    }
}
</script>

<template>
    <div class="container-fluid">
        <div class="row">
            <div class="col-lg-8">
                <div id="viewport" style="width: 100%; height: 500px;"></div>
                <div class="amino-acid-sequence" style="width: 100%; color: black; background-color: #f8f8f8; padding: 10px;">
                    <span v-for="(acid, index) in sequence" :key="index"
                          :class="{ 'selected': selectedResidues.has(index) }"
                          @click="toggleResidueSelection(index)">
                        {{ acid }}
                    </span>
                </div>
            </div>
            <div class="col-lg-4">
                <h4>Options</h4>
                <select class="form-select form-select-md mb-3" aria-label="Select a representation" @change="setView($event)">
                    <option v-for="viewKey in viewsItems" :key="viewKey" :value="viewKey">
                        {{ views[viewKey].title }}
                    </option>
                </select>
                <button type="button" @click="downloadPdbFile()"
                    class="btn btn-primary padding-top-button-group download-pdb-button">Download .pdb
                </button>
            </div>
        </div>
        <div class="row">
            <div class="col-12" style="margin-top: 10px;">
                <button class="btn btn-danger" @click="$emit('close')">Close</button>
            </div>
        </div>
    </div>
</template>

<style>
.amino-acid-sequence {
    overflow-x: auto;
    white-space: nowrap;
    font-family: monospace;
    margin-bottom: 10px;
    border: 1px solid #ddd; /* Optional: adds a border around the sequence */
    border-radius: 5px; /* Optional: rounds the corners of the border */
}

.amino-acid-sequence span {
    display: inline-block;
    padding: 2px 5px;
    cursor: pointer;
    border-radius: 3px;
}

.amino-acid-sequence span.selected {
    background-color: blue;
    color: white;
}

.loading-indicator {
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 10px;
}
</style>
