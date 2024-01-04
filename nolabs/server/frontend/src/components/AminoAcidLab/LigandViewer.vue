<script>
export default {
    props: ['api', 'experiment', 'ligand'],
    data() {
        const views = {
            // Define different views suitable for ligands
            ballsAndSticks: { key: 'ball+stick', title: 'Balls and Sticks' },
            spacefill: { key: 'spacefill', title: 'Spacefill' },
            // ... other views if needed ...
        };
        return {
            views: views,
            viewsItems: Object.keys(views),
            selectedView: views.ballsAndSticks, // Default view
            stage: null,
            sdfComponent: null,
            isLoading: false,
        };
    },
    methods: {
        setView(evt) {
            const viewKey = evt.target.value;
            this.selectedView = this.views[viewKey];
            this.sdfComponent.removeAllRepresentations();
            this.sdfComponent.addRepresentation(this.selectedView.key);
        },
        render() {
            setTimeout(() => {
                document.getElementById('viewport').innerHTML = '';
                this.stage = new NGL.Stage("viewport");
                this.stage.setParameters({ backgroundColor: 'white' });
                const sdfContentBlob = new Blob([this.ligand.sdf], { type: 'chemical/x-mdl-sdfile' });
                const sdfFile = new File([sdfContentBlob], 'ligand.sdf', { type: 'chemical/x-mdl-sdfile' });
                this.stage.loadFile(sdfFile, { defaultRepresentation: true }).then((component) => {
                    this.sdfComponent = component;
                });
            }, 200);
        },
        downloadSdfFile() {
            const filename = this.ligand.metadata.name + '.sdf';
            const blob = new Blob([this.ligand.sdf], { type: 'chemical/x-mdl-sdfile' });
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
    }
}
</script>

<template>
    <div class="container-fluid">
        <div class="row">
            <div class="col-lg-8">
                <div id="viewport" style="width: 100%; height: 500px;"></div>
            </div>
            <div class="col-lg-4">
                <h4>Options</h4>
                <select class="form-select form-select-md mb-3" aria-label="Select a representation" @change="setView($event)">
                    <option v-for="viewKey in viewsItems" :key="viewKey" :value="viewKey">
                        {{ views[viewKey].title }}
                    </option>
                </select>
                <button type="button" @click="downloadSdfFile"
                    class="btn btn-primary padding-top-button-group">Download .sdf
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
