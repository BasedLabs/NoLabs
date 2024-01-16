<template>
    <q-page-container>
        <q-drawer show-if-above v-model="drawer" side="left" bordered>
            <q-list bordered class="rounded-borders">
                <q-item-label header>Targets</q-item-label>
                <q-expansion-item v-for="target in targets" :key="target.target_id" :label="target.target_name"
                    class="q-ma-none">
                    <q-list>
                        <q-item v-for="ligand in getLigandsForTarget(target.target_id)" :key="ligand.ligand_id" clickable>
                            {{ ligand.ligand_name }}
                        </q-item>
                    </q-list>
                </q-expansion-item>
            </q-list>
        </q-drawer>


        <q-page padding>
            <div v-if="loading">
                <q-spinner color="primary" />
            </div>
            <!-- Main content area, display experiment details or other information here -->
        </q-page>
    </q-page-container>
</template>
<script lang="ts">
import { defineComponent, computed, ref } from 'vue';
import { useDrugDiscoveryStore } from './storage';
import { useRoute } from 'vue-router';
import { QSpinner, QDrawer, QList, QItem, QItemLabel, QExpansionItem } from 'quasar';

export default defineComponent({
    name: 'ExperimentDetail',
    components: {
        QSpinner,
        QDrawer,
        QList,
        QItem,
        QItemLabel,
        QExpansionItem,
    },
    setup() {
        const store = useDrugDiscoveryStore();
        const route = useRoute();
        const experimentId = route.params.experimentId as string;
        const loading = ref(true);
        const drawer = ref(true);

        store.fetchTargetsForExperiment(experimentId).then(() => {
            loading.value = false;
        }).catch((error) => {
            console.error('Error fetching targets:', error);
            loading.value = false;
        });

        const getLigandsForTarget = (targetId: string) => {
            store.fetchLigandsForTarget(experimentId, targetId);
            return store.ligands;
        };

        const targets = computed(() => store.targets);

        return {
            targets,
            getLigandsForTarget,
            loading,
            drawer,
        };
    },
});
</script>
<style>
.rounded-borders {
    border-radius: 8px;
}
</style>
  