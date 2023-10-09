<script setup>

import { ref, reactive, onMounted } from 'vue';
import store from '../../storage';

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

const render = (experiment, bindingIndex) => {
    if(bindingIndex >= experiment.data.length){
        console.log('No experiment data');
        return;
    }
    const drugTargetData = experiment.data[bindingIndex];
    
    setTimeout(() => {
        // Render 3d structure
        document.getElementById('viewport').innerHTML = '';
        cleanStage();
        stage = new NGL.Stage("viewport");
        stage.setParameters({ backgroundColor: 'white' });
        const proteinFileContentBlob = new Blob([drugTargetData.pdb], { type: 'text/plain' });
        const ligandFileContentBlob = new Blob([drugTargetData.sdf], { type: 'text/plain' });
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

const renderTable = (experiment) => {
    const experimentData = experiment.data;
    // Render table
    const tableData = [];
    for (let [index, val] of experimentData.entries()) {
        tableData.push({ id: index, ligand: val.ligandName, protein: val.proteinName, affinity: val.affinity });
    }

    const rowSelect = (rowData) => {
        render(experiment, rowData.id);
    }

    $('#ligandProteinTable').bootstrapTable({
        data: tableData,
        onClickRow: (row, el, field) => {
            $('#ligandProteinTable tr').removeClass('active');
            $(el).addClass('active');
            rowSelect(row)
        }
    });

    $("#ligandProteinTableSearch").on("keyup", function () {
        const value = $(this).val().toLowerCase();
        $("#ligandProteinTable tr").filter(function () {
            $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
        });
    });
}

onMounted(() => {
    renderTable(store.state.drugTargetData.experiment);
    render(store.drugTargetData.experiment, 0);
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
        <div class="row">
            <div class="col-md-12 container">
                <h4>Ligand-protein pairs</h4>
                <label for="ligandProteinTableSearch">Search for ligand-protein pair, click to render:</label><input
                    class="form-control" id="ligandProteinTableSearch" type="text" placeholder="Search..">
                <br>
                <table style="max-height: 200px; overflow-y:scroll;" class="table table-bordered table-striped"
                    id="ligandProteinTable">
                    <thead>
                        <tr>
                            <th data-field="ligand">Ligand</th>
                            <th data-field="protein">Protein</th>
                            <th data-field="affinity">Affinity</th>
                        </tr>
                    </thead>
                </table>
            </div>
            <div class="col-md-2"></div>
        </div>
    </div>
</template>