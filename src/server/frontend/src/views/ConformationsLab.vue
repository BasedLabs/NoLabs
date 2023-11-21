<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import DrugTargetLabFormView from './DrugTargetLabFormView.vue';
import DrugTargetLabInferenceView from './DrugTargetLabInferenceView.vue';

export default {
    data() {
        return {
            experimentsState: store.state.drugTarget,
            experimentsApi: api.drugTarget
        }
    },
    components: {
        ExperimentsView,
        DrugTargetLabFormView,
        DrugTargetLabInferenceView
    }
}
</script>

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
            experimentsState: store.state.conformations,
            experimentsApi: api.conformations,
            views: views,
            viewsItems: Object.keys(views),
            selectedView: views.default,
            stage: null,
            pdbComponent: null,
            ligandComponent: null
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
                const experiment = store.state.conformations.experiment;
                document.getElementById('viewport').innerHTML = '';
                this.cleanStage();
                const stage = new NGL.Stage("viewport");
                stage.setParameters({ backgroundColor: 'black' });
                const proteinFileContentBlob = new Blob([experiment.data], { type: 'text/plain' });
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
    },
    components: {
        ExperimentsView
    }
}
</script>

<template>
    <ExperimentsView :state="experimentsState" :api="experimentsApi">
        <template v-slot:labTitle>
            <h4>Drug target lab</h4>
        </template>
        <template v-slot:labForm="labForm">
            <DrugTargetLabFormView :onFormSubmit="labForm.onFormSubmit"/>
        </template>
        <template v-slot:lab="lab">
            <DrugTargetLabInferenceView :experiment="lab.experiment"/>
        </template>
    </ExperimentsView >
</template>