<script>
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
            pdbContent: null,
            views: views,
            viewsItems: Object.keys(views),
            selectedView: views.default,
            stage: null,
            pdbComponent: null,
            ligandComponent: null
        }
    },
    methods: {
        loadPdb(event) {
            const reader = new FileReader();
            reader.addEventListener(
                "load",
                () => {
                    this.pdbContent = reader.result;
                    this.render()
                },
                false,
            );
            reader.readAsText(event.target.files[0]);
        },
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
                this.stage.setParameters({ backgroundColor: 'black' });
                const proteinFileContentBlob = new Blob([this.pdbContent], { type: 'text/plain' });
                const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
                this.stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
                    this.pdbComponent = component;
                });
            }, 200);
        }
    }
}
</script>

<template>
    <div class="container-fluid">
        <div class="row height-100">
            <div class="col-md-3 text-center justify-content-center">
                <h4 v-if="pdbContent != null">Options</h4>
                <select v-if="pdbContent != null" class="form-select form-select-md" aria-label="Select a representation"
                    @change="setView($event)">
                    <option v-for="viewKey in viewsItems" :key="viewKey" :value="viewKey">
                        {{ views[viewKey].title }}
                    </option>
                </select>
                <input type="file" ref="file" style="display: none" @change="loadPdb" accept=".pdb" />
                <button type="button" @click="$refs.file.click()"
                    class="btn btn-primary padding-top-button-group download-pdb-button">Upload .pdb
                </button>
            </div>
            <div class="col-md-9 text-center">
                <h4>Pdb viewer</h4>
                <div id="viewport" style="width: 100%; min-height: 80vh;"></div>
            </div>
        </div>

    </div>
</template>