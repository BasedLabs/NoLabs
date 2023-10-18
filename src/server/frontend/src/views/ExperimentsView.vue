<script>


export default {
    props: ['experimentsState', 'experimentsApi'],
    data() {
        return {
            experimentsState: this.experimentsState,
            experimentsApi: this.experimentsApi
        }
    },
    methods: {
        loadExperiments: async () => {
            await this.api.getAllExperiments();
        },
        async selectExperiment(experiment) {
            this.$loading.show();
            await this.api.loadExperiment(experiment);
            this.$loading.hide();
        },
        addExperiment: async () => {
            this.api.addExperiment();
        },
        deleteExperiment: async (experiment) => {
            await this.api.deleteExperiment(experiment);
            await loadExperiments();
        },
        async onFormSubmit(data) {
            await this.api.inference(data.target);
            await this.selectExperiment()
        },
        isCurrentExperiment(experiment){
            return this.state.experiment && (this.state.experiment.name == experiment.name)
        }
    },
    async mounted() {
        await this.loadExperiments();
    },
    computed: {
        exprimentNotEmpty() {
            return !!this.state.experiment && this.state.experiment.data && this.state.experiment.data.length > 0;
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
                <li v-for="experiment in drugTarget.experiments"
                    class="list-group-item list-group-item-action d-flex justify-content-between align-items-center"
                    :class="isCurrentExperiment(experiment) ? 'active' : ''"
                    @click.stop="selectExperiment(experiment)">
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
                        <slot name="labForm"></slot>
                    </div>
                </div>
                <div id="resultContainer" v-if="exprimentNotEmpty">
                    <div class="tab-content" id="nav-tabContent">
                        <div class="tab-pane fade mt-1 show active" role="tabpanel">
                            <slot name="lab" :key="state.experiment"></slot>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>
