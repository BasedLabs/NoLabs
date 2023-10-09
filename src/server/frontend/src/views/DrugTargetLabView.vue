<script setup>
import DrugTarget from '../components/DrugTargetDiscoveryLab/DrugTarget.vue';
import { onMounted, reactive, defineEmits, ref } from 'vue';
import store, {api} from '../storage';

const experiments = store.state.drugTargetData.experiments;
const experiment = store.state.drugTargetData.experiment;

const experimentLoaded = Object.keys(experiment).length > 0;

const loadExperiments = async () => {
    await api.drugTargetDiscovery.getAllExperiments();
    console.log(experiment.length);
}

const selectExperiment = async (experiment) => {
    await api.drugTargetDiscovery.loadExperiment(experiment.id);
}

const addExperiment = async () => {
    await api.drugTargetDiscovery.addExperiment();
    await loadExperiments();
}

const deleteExperiment = async (name) => {
    await api.drugTargetDiscovery.deleteExperiment();
    await loadExperiments();
}

onMounted(async () => {
    await loadExperiments();
})
</script>

<template>
    <div class="row gy-1">
        <div class="col-md-2 experiments-col">
            <ul class="list-group">
                <li class="list-group-item d-flex justify-content-between align-items-center">
                    <button type="button" @click="addExperiment()"
                        class="btn btn-outline-success add-experiments-button">Add</button>
                </li>
                <li v-for="experiment in experiments"
                    class="list-group-item list-group-item-action d-flex justify-content-between align-items-center"
                    @click="selectExperiment(experiment)"
                    >
                    {{ experiment.name }}
                    <span class="badge bg-danger rounded-pill btn btn-outline-danger btn-sm btn-link text-decoration-none"
                        @click="deleteExperiment(experiment.name)">X</span>
                </li>
            </ul>
        </div>
        <div class="col-md-10">
            <div class="text-center m-5">
                <h4>Drug target discovery lab</h4>
                <div class="row">
                    <div class="col-md-12">
                        <form enctype="multipart/form-data" id="inferenceInputForm">
                            <div class="row justify-content-center">
                                <div class="col-md-6">
                                    <label for="sdfFileInput" class="col-form-label fs-4">Molecules (.sdf format)</label>
                                    <input type="file" multiple="multiple" accept="text/x-smi" class="form-control"
                                        id="sdfFileInput" name="smilesFileInput">
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

                <div class="text-center invisible" id="spinner">
                    <p class="inference-components-title">Processing. It can take more than 10 mins depending on your
                        PC</p>
                    <div class="spinner-border" role="status">
                        <span class="visually-hidden">Processing</span>
                    </div>
                </div>
                <div id="resultContainer" v-if="experimentLoaded">
                    <div class="tab-content" id="nav-tabContent">
                        <div class="tab-pane fade mt-1 show active" role="tabpanel">
                            <DrugTarget />
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>
