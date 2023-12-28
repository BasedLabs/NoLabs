<script>

export default {
    props: ['experiment', 'api'],
    components: {
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
    </div>
</template>
