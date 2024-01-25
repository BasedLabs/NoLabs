<template>
    <div>
        <q-page-container>
            <q-toolbar>
                <!-- Tab Menu -->
                <q-tabs aligned="justify" class="text-grey bg-white full-width" active-color="primary"
                    indicator-color="primary">
                    <q-tab name="tab1" label="Experiment setup" @click="openExperimentSetup" />
                    <q-tab name="tab2" label="Running jobs" @click="openRunningJobs" />
                    <q-tab name="tab3" label="Completed" />
                    <!-- Add more tabs as needed -->
                </q-tabs>
            </q-toolbar>
            <router-view></router-view> <!-- This tag renders the active route component -->
        </q-page-container>
    </div>
</template>
  
<script>
export default {
    name: 'ExperimentNavigation',
    data() {
        return {
            experimentId: null,
        };
    },
    methods: {
        openExperimentSetup() {
            this.$router.push({ name: 'ExperimentSetup', params: { expId: this.experimentId } });
        },
        openRunningJobs() {
            this.$router.push({ name: 'RunningJobs', params: { expId: this.experimentId } });
        }
    },
    mounted() {
        this.experimentId = this.$route.params.experimentId;
    },
};
</script>