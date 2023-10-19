<script>
export default {
    props: ['state', 'api'],
    methods: {
        async loadExperiments () {
            await this.api.getAllExperiments();
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
            await this.api.deleteExperiment(experiment);
            await loadExperiments();
        },
        async onFormSubmit(data) {
            await this.api.inference({form: data.target, experiment: this.state.experiment});
            await this.selectExperiment()
        },
        isCurrentExperiment(experiment) {
            return this.state.experiment && (this.state.experiment.name == experiment.name)
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
                <li v-for="experiment in state.experiments"
                    class="list-group-item list-group-item-action d-flex justify-content-between align-items-center"
                    :class="isCurrentExperiment(experiment) ? 'active' : ''" @click.stop="selectExperiment(experiment)">
                    {{ experiment.name }}
                    <span class="badge bg-danger rounded-pill btn btn-outline-danger btn-sm btn-link text-decoration-none"
                        @click.stop="deleteExperiment(experiment)">X</span>
                </li>
            </ul>
        </div>
        <div class="col-md-10">
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
