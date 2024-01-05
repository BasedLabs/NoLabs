<script>
import DrugTargetResultsList from './DrugTargetResultsList.vue';


export default {
    props: ['experiment', 'api'],
    data() {
        return {
            experimentRunning: false
        };
    },
    components: {
        DrugTargetResultsList,
    },
    methods: {
        runExperiment() {
            this.api.inference(this.experiment);
            this.experimentRunning = true;
        },
    },
    mounted() {
    }
}
</script>

<template>
    <div>
        <button class="btn" 
                :class="experimentRunning ? 'btn-secondary' : 'btn-success'" 
                @click="runExperiment" 
                :disabled="experimentRunning">
            Run Experiment
        </button>
        <DrugTargetResultsList :experiment="experiment" :api="api"></DrugTargetResultsList>
    </div>
</template>
