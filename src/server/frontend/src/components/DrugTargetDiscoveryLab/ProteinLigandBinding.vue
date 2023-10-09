<script setup>
import store from '../../storage.js';

import { ref, reactive, onMounted } from 'vue';

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
const viewsItems = Object.keys(views);

let stage;
let pdbComponent;
let ligandComponent;
let selectedView = ref(views.default);

const setView = (viewKey) => {
    selectedView = views[viewKey];
    if (selectedView === views.default) {
        render();
        return;
    }
    pdbComponent.removeAllRepresentations();
    pdbComponent.addRepresentation(selectedView.key);
    ligandComponent.removeAllRepresentations();
    ligandComponent.addRepresentation(selectedView.key);
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
        const proteinFileContentBlob = new Blob([drugTargetData.inference.pdb], { type: 'text/plain' });
        const ligandFileContentBlob = new Blob([drugTargetData.inference.sdf], { type: 'text/plain' });
        const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
        const sdfFile = new File([ligandFileContentBlob], 'ligand.sdf', { type: 'text/plain' });
        stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
            pdbComponent = component;
        });

        stage.loadFile(sdfFile, { defaultRepresentation: true }).then((component) => {
            ligandComponent = component;
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
    <div class="row">
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
        </div>
    </div>
</template>