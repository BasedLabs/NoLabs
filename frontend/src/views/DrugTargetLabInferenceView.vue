<script>
import store, { api } from '../storage';

export default {
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
            component: null,
            sdfComponent: null,
            selectedRowId: 0
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
            this.component.removeAllRepresentations();
            this.component.addRepresentation(this.selectedView.key);
        },
        cleanStage() {
            if (this.stage)
                this.stage.removeAllComponents();
        },
        render() {
            const drugTarget = this.experiment.data;

            setTimeout(() => {
                // Render 3d structure
                document.getElementById('viewport').innerHTML = '';
                this.cleanStage();
                this.stage = new NGL.Stage("viewport");
                this.stage.setParameters({ backgroundColor: 'white' });
                const proteinFileContentBlob = new Blob([drugTarget.pdb], { type: 'text/plain' });
                const ligandFileContentBlob = new Blob([drugTarget.sdf], { type: 'text/plain' });
                const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
                const sdfFile = new File([ligandFileContentBlob], 'ligand.sdf', { type: 'text/plain' });
                this.stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
                    this.component = component;
                    this.stage.loadFile(sdfFile, { defaultRepresentation: true }).then((sdfComponent) => {
                        this.sdfComponent = sdfComponent;
                    });
                });
            }, 200);
        },
        zoomInLigand() {
            this.sdfComponent.autoView();
        },
        zoomInProtein() {
            this.component.autoView();
        },
        async downloadPdbFile(evt) {
            const response = await api.drugTarget.downloadCombinedPdb(this.experiment, this.experiment.metaData['selectedProteinId'], this.experiment.metaData['selectedLigandId']);
            const filename = this.experiment.metaData['selectedProteinId'] + this.experiment.metaData['selectedLigandId'] + '_docking.pdb';
            const blob = new Blob([response.data.pdb], { type: 'text/plain' });
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
    mounted() {
        this.render();
    },
    computed: {
        experiment() {
            return store.state.drugTarget.experiment;
        }
    }
}
</script>

<template>
    <div class="row">
        <div class="col-md-8">
            <div id="viewport" style="width: 100%; height: 500px;"></div>
        </div>
        <div class="col-md-4" style="margin-top: 10px;">
            <h4>Options</h4>
            <select class="form-select form-select-md mb-3" aria-label="Select a representation" @change="setView($event)">
                <option v-for="viewKey in viewsItems" :key="viewKey" :value="viewKey">
                    {{ views[viewKey].title }}
                </option>
            </select>
            <button type="button" @click="zoomInLigand()" style="width: 100%"
                class="btn btn-primary padding-top-button-group">Zoom in ligand
            </button>
            <button type="button" @click="zoomInProtein()" style="width: 100%"
                class="btn btn-primary padding-top-button-group">Zoom in protein
            </button>
            <button type="button" @click="downloadPdbFile()" style="width: 100%"
                class="btn btn-primary padding-top-button-group">Download structure .pdb
            </button>
        </div>
    </div>
</template>