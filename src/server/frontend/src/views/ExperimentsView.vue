<script>
export default {
    props: ['state', 'api'],
    methods: {
        async loadExperiments () {
            const loader = this.$loading.show();
            await this.api.getAllExperiments();
            loader.hide();
        },
        async selectExperiment(experiment) {
            const loader = this.$loading.show();
            await this.api.loadExperiment(experiment);
            loader.hide();
        },
        async addExperiment() {
            this.api.addExperiment();
        },
        async deleteExperiment(experiment) {
            const loader = this.$loading.show();
            await this.api.deleteExperiment(experiment);
            loader.hide();
            await loadExperiments();
        },
        async onFormSubmit(data) {
            const loader = this.$loading.show();
            await this.api.inference({form: data.target, experiment: this.state.experiment});
            loader.hide();
        },
        isCurrentExperiment(experiment) {
            return this.state.experiment && (this.state.experiment.name == experiment.name)
        },
        showTooltip(evt, text) {
            let tooltip = document.getElementById("tooltip");
            tooltip.innerHTML = text;
            tooltip.style.display = "block";
            tooltip.style.left = evt.pageX + 10 + 'px';
            tooltip.style.top = evt.pageY + 10 + 'px';
        },
        hideTooltip() {
            var tooltip = document.getElementById("tooltip");
            tooltip.style.display = "none";
        },
        async changeExperimentName(evt, experiment) {
            experiment.name = evt.target.value;
            await this.api.changeExperimentName(experiment);
        }
    },
    async mounted() {
        await this.loadExperiments();
    },
    computed: {
        experimentEmpty() {
            if(!this.state.experiment || !this.state.experiment.data)
                return true;

            if(typeof this.state.experiment.data === 'Array')
                return this.state.experiment.data.length === 0;

            if(Object.keys(this.state.experiment.data).length > 0)
            {
                console.log('NOT EMPTY');
            }

            return Object.keys(this.state.experiment.data).length === 0;
        },
        experimentSelected() {
            return !!this.state.experiment;
        }
    }
}
</script>

<template>
    <div class="row gy-1">
        <div class="col-md-2 experiments-col">
            <ul class="list-group">
                <li class="list-group-item d-flex justify-content-between align-items-center">
                    <button type="button" @click.stop="addExperiment()"
                        class="btn btn-outline-success add-experiments-button">Add</button>
                </li>
                <div id="tooltip" display="none" style="position: absolute; display: none; z-index: 10"></div>
                <li v-for="experiment in state.experiments"
                    class="list-group-item list-group-item-action d-flex justify-content-between align-items-center"
                    :class="isCurrentExperiment(experiment) ? 'active' : ''">
                    <input class="form-control" :value="experiment.name" @input="changeExperimentName($event, experiment)"/>
                    <i class="bi bi-download btn btn-sm btn-outline-success" type="button" 
                        @click.stop="selectExperiment(experiment)"
                        @mouseover="showTooltip($event, 'Load into solution')"
                        @mouseout="hideTooltip()"></i>
                    <i class="bi bi-x-circle btn btn-sm btn-outline-danger"
                        @click.stop="deleteExperiment(experiment)"
                        @mouseover="showTooltip($event, 'Delete from solution')"
                        @mouseout="hideTooltip()"></i>
                </li>
            </ul>
        </div>
        <div class="col-md-8">
            <div class="text-center m-5">
                <slot name="labTitle"></slot>
                <div class="row" v-if="experimentSelected">
                    <div class="col-md-12">
                        <slot name="labForm" :onFormSubmit="onFormSubmit"></slot>
                    </div>
                </div>
                <div class="row" v-else>
                    <div class="col-md-12">
                        <hr/>
                        <h4>Add or select experiment</h4>
                    </div>
                </div>
                <div id="resultContainer" v-if="!experimentEmpty">
                    <div class="tab-content" id="nav-tabContent">
                        <div class="tab-pane fade mt-1 show active" role="tabpanel">
                            <slot name="lab" :key="state.experiment.id" :experiment="state.experiment"></slot>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>
