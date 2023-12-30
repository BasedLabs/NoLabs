<script>
import DrugTargetResultsList from './DrugTargetResultsList.vue';


export default {
    props: ['experiment', 'api'],
    components: {
        DrugTargetResultsList,
    },
    data() {
        return {
        };
    },
    methods: {
        runExperiment() {
            this.api.inference(this.experiment);
        },
    },
    mounted() {
    }
}
</script>

<template>
    <div>
        <button class="btn btn-success" @click="runExperiment"> Run Experiment </button>
        <DrugTargetResultsList :experiment="experiment" :api="api"></DrugTargetResultsList>
    </div>
</template>
