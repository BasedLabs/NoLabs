<script>
import DrugTarget from '../components/DrugTargetDiscoveryLab/DrugTarget.vue';
import store, { api } from '../storage';

export default {
    data() {
        return {
            drugTargetData: store.state.drugTargetData,
        }
    },
    methods: {
        loadExperiments: async () => {
            await api.drugTargetDiscovery.getAllExperiments();
        },
        async selectExperiment(experiment) {
            await api.drugTargetDiscovery.loadExperiment(experiment.name);
        },
        addExperiment: async () => {
            api.drugTargetDiscovery.addExperiment();
        },
        deleteExperiment: async (experiment) => {
            await api.drugTargetDiscovery.deleteExperiment(experiment);
            await loadExperiments();
        },
        async onFormSubmit(data) {
            await api.drugTargetDiscovery.inference(data.target);
            await this.selectExperiment()
        },
        isCurrentExperiment(experiment){
            return this.drugTargetData.experiment && (this.drugTargetData.experiment.name == experiment.name)
        }
    },
    async mounted() {
        await this.loadExperiments();
    },
    computed: {
        exprimentNotEmpty() {
            return !!this.drugTargetData.experiment && this.drugTargetData.experiment.data && this.drugTargetData.experiment.data.length > 0;
        },
        experimentSelected() {
            return !!this.drugTargetData.experiment;
        }
    },
    components: {
        DrugTarget
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
                <li v-for="experiment in drugTargetData.experiments"
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
                <h4>Drug discovery lab</h4>
                <div class="row">
                    <div class="col-md-12">
                        <form enctype="multipart/form-data" id="inferenceInputFormDrugDiscovery"
                            v-on:submit.prevent="onFormSubmit" v-if="experimentSelected">
                            <div class="row justify-content-center">
                                <div class="col-md-6">
                                    <label for="sdfFileInput" class="col-form-label fs-4">Molecules (.sdf format)</label>
                                    <input type="file" multiple="multiple" accept="text/x-smi" class="form-control"
                                        id="sdfFileInput" name="sdfFileInput">
                                </div>

                                <div class="col-md-6">
                                    <label for="proteinFileInput" class="col-form-label fs-4">Proteins (.pdb format)</label>
                                    <input type="file" accept="text/x-pdb" class="form-control" id="proteinFileInput"
                                        name="proteinFileInput">
                                </div>
                            </div>
                            <div class="row justify-content-center m-4">
                                <div class="col-md-4">
                                    <div class="btn-group mt-4" role="group"
                                        aria-label="Basic checkbox toggle button group">
                                        <button type="submit" id="submitInference" class="btn btn-lg btn-primary">
                                            Submit
                                        </button>
                                    </div>
                                </div>
                            </div>
                        </form>
                    </div>
                </div>
                <div id="resultContainer" v-if="exprimentNotEmpty">
                    <div class="tab-content" id="nav-tabContent">
                        <div class="tab-pane fade mt-1 show active" role="tabpanel">
                            <DrugTarget :key="drugTargetData.experiment.name" />
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>
