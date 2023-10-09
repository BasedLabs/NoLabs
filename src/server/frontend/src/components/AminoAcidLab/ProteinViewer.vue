<script setup>
import store from '../../storage.js';

import { ref, onMounted } from 'vue'

const downloadPdbFile = () => {
    const filename = 'protein.pdb';
    const blob = new Blob([store.state.aminoAcidData.inference.folding], { type: 'text/plain' });
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

const views = {
    default: {key: 'default', title: 'Default representation'},
    cartoon: {key: 'cartoon', title: 'Cartoon'},
    backbone: {key: 'backbone', title: 'Backbone'},
    ballsAndSticks: {key: 'ball+stick', title: 'Balls and sticks'},
    contact: {key: 'contact', title: 'Contact'},
    helixorient: {key: 'helixorient', title: 'Helixorient'},
    hyperball: {key: 'hyperball', title: 'Hyperball'},
    licorice: {key: 'licorice', title: 'Licorice'},
    ribbon: {key: 'ribbon', title: 'Ribbon'},
    rope: {key: 'rope', title: 'Rope'},
    surface: {key: 'surface', title: 'Surface'},
    spacefill: {key: 'spacefill', title: 'Spacefill'},
    unitcell: {key: 'unitcell', title: 'Unitcell'}
};
const viewsItems = Object.keys(views);

let nglComponent;
let stage;
let selectedView = ref(views.default);

const setView = (viewKey) => {
    selectedView = views[viewKey];
    if (selectedView === views.default) {
        render();
        return;
    }
    nglComponent.removeAllRepresentations();
    nglComponent.addRepresentation(selectedView.key);
}

const cleanStage = () => {
    if (stage)
        stage.removeAllComponents();
}

const render = () => {
    setTimeout(() => {
        document.getElementById('viewport').innerHTML = '';
        cleanStage();
        stage = new NGL.Stage("viewport");
        stage.setParameters({ backgroundColor: 'white' });
        const proteinFileContentBlob = new Blob([store.state.aminoAcidData.inference.folding], { type: 'text/plain' });
        const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
        stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
            nglComponent = component;
        });
    }, 500);

}

onMounted(() => {
    render();
});

defineExpose({
    render
});

</script>

<template>
    <div class="row mt-2">
        <div class="col-md-8">
            <div id="viewport" style="width: 100%; height: 500px;"></div>
        </div>
        <div class="col-md-4">
            <h4>Options</h4>
            <select class="form-select form-select-md mb-3" aria-label="Select a representation">
                <option v-for="viewKey in viewsItems" @click="setView(viewKey)" :selected="selectedView === view">
                    {{ views[viewKey].title }}
                </option>
            </select>
            <button type="button" @click="downloadPdbFile()" class="btn btn-primary padding-top-button-group download-pdb-button">Download .pdb
            </button>
        </div>
    </div>
</template>
