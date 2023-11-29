<script>
export default {
    props: ['experiment'],
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
            pdbComponent: null
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
                const experiment = this.experiment;
                document.getElementById('viewport').innerHTML = '';
                this.cleanStage();
                const stage = new NGL.Stage("viewport");
                stage.setParameters({ backgroundColor: 'white' });
                const proteinFileContentBlob = new Blob([experiment.data.pdb], { type: 'text/plain' });
                const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
                stage.loadFile(proteinFile, { defaultRepresentation: true, asTrajectory: true }).then(function(component) {
                    var traj = component.addTrajectory().trajectory
                    var player = new NGL.TrajectoryPlayer(traj, {
                        step: 1,
                        timeout: 700,
                        interpolateStep: 100,
                        start: 0,
                        end: traj.numframes,
                        interpolateType: "linear",
                        mode: "loop",
                        direction: "bounce"
                    });
                    player.play();

                    component.autoView();
                });
            }, 200);
        }
    },
    mounted() {
        this.render();
    }
}
</script>

<template>
    <div class="row mt-2">
        <div class="col-md-8">
            <div id="viewport" style="width: 100%; height: 500px;"></div>
        </div>
        <div class="col-md-4" style="align-self: center; margin-top: 10px;">
            <h4>Options</h4>
            <select class="form-select form-select-md mb-3" aria-label="Select a representation" @change="setView($event)">
                <option v-for="viewKey in viewsItems" :key="viewKey" :value="viewKey">
                    {{ views[viewKey].title }}
                </option>
            </select>
        </div>
    </div>
</template>
