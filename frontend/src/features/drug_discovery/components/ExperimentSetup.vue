<template>
    <q-page>
        <!-- Main Page Content -->
        <q-page>
            <div v-if="loading">
                <q-spinner color="primary" />
            </div>
            <div v-else>
                <TargetsList :experimentId="this.experimentId" :originalTargets="this.targets"/>
            </div>
        </q-page>
    </q-page>
</template>

<script>
import { useDrugDiscoveryStore } from '../storage';
import { QSpinner } from 'quasar';
import { useRoute } from 'vue-router';
import TargetsList from './targets/TargetsList.vue'

export default {
    name: 'ExperimentSetup',
    components: {
        TargetsList
    },
    data() {
        return {
            loading: true,
            targets: null
        };
    },
    mounted() {
        const store = useDrugDiscoveryStore();
        const route = useRoute();
        this.experimentId = route.params.experimentId;
        store.fetchTargetsForExperiment(this.experimentId).then(() => {
            this.targets = store.targets;
            this.loading = false;
        }).catch((error) => {
            console.error('Error fetching targets:', error);
            this.loading = false;
        });
    }
};
</script>
