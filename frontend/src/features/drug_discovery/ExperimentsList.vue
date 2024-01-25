<template>
    <q-page padding>
        <div class="q-pa-md text-center">
            <q-btn color="green" label="Add Experiment" @click="addNewExperiment" />
        </div>

        <div v-if="loading">
            <q-spinner color="primary" />
        </div>
        <div v-else class="experiment-list">
            <div v-for="experiment in experiments" :key="experiment.experiment_id" class="experiment-tab row q-mb-md">
                <div @click="openExperiment(experiment.experiment_id)">
                    <div class="col experiment-info" @click="openExperiment(experiment.experiment_id)">
                        <div class="experiment-name">{{ experiment.experiment_name }}</div>
                        <div class="experiment-date">{{ formatDate(experiment.experiment_date) }}</div>
                    </div>
                    <div class="col-auto">
                        <q-btn flat icon="more_vert" @click.stop="showMenu(experiment.experiment_id)">
                            <q-menu>
                                <q-list>
                                    <q-item clickable @click="changeExperimentName(experiment.experiment_id)">
                                        <q-item-section>Change Name</q-item-section>
                                    </q-item>
                                    <q-item clickable @click="deleteExperiment(experiment.experiment_id)">
                                        <q-item-section>Delete</q-item-section>
                                    </q-item>
                                </q-list>
                            </q-menu>
                        </q-btn>
                    </div>
                </div>
            </div>
        </div>
    </q-page>
</template>
  
<script lang="ts">
import { defineComponent, ref, computed } from 'vue';
import { useRouter } from 'vue-router';
import { useDrugDiscoveryStore } from './storage';
import { QSpinner, QBtn, QMenu, QList, QItem, QItemSection, useQuasar } from 'quasar';

export default defineComponent({
    name: 'ExperimentsList',
    components: {
        QSpinner,
        QBtn,
        QMenu,
        QList,
        QItem,
        QItemSection,
    },
    setup() {
        const store = useDrugDiscoveryStore();
        const loading = ref(true);
        const router = useRouter();
        const $q = useQuasar();

        store.fetchExperiments().then(() => {
            loading.value = false;
        }).catch((error) => {
            console.error('Error fetching experiments: ', error);
            loading.value = false;
        });

        const experiments = computed(() => store.experiments);

        const addNewExperiment = async () => {
            try {
                const newExperiment = await store.createExperiment();
                $q.notify(`Successfully added a new experiment`);
            } catch (error) {
                console.error('Error adding experiment:', error);
                $q.notify('Error adding experiment');
            }
        };


        const deleteExperiment = (experimentId: string) => {
            store.removeExperiment(experimentId);
        };

        const showMenu = (experiment_id: string) => {
            // Logic to handle menu display
        };

        const changeExperimentName = (experiment_id: string) => {
            // Logic to change the experiment's name
        };


        const loadMore = () => {
            // Placeholder for infinite scroll logic
            // Since all experiments are loaded at once, no additional loading is needed
        };

        const openExperiment = (experimentId: string) => {
            router.push({ name: 'ExperimentNavigation', params: { experimentId } });
        };

        const formatDate = (dateString: string) => {
            if (!dateString) return '';
            const date = new Date(dateString);
            return date.toLocaleDateString();
        };

        return {
            experiments,
            loading,
            addNewExperiment,
            showMenu,
            changeExperimentName,
            deleteExperiment,
            formatDate,
            openExperiment,
        };
    },
});
</script>
  
<style>
.experiment-list {
    max-height: 500px;
    /* Adjust as needed */
    overflow-y: auto;
}

.experiment-tab {
    border: 1px solid #E0E0E0;
    border-radius: 8px;
    padding: 15px;
    display: flex;
    justify-content: space-between;
    align-items: center;
}

.experiment-info {
    flex-grow: 1;
    cursor: pointer;
}

.experiment-name {
    font-weight: bold;
    font-size: 1.1em;
}

.experiment-date {
    color: #757575;
}
</style>