<template>
    <q-page>
        <!-- Main Page Content -->
        <q-page>
            <div v-if="loading">
                <q-spinner color="primary" />
            </div>
            <div class="q-pa-md">
                <q-table title="Queue" :rows="running_jobs" :columns="columns" row-key="id" v-model:expanded="expanded">

                    <template v-slot:header="props">
                        <q-tr :props="props">
                            <q-th auto-width />
                            <q-th v-for="col in props.cols" :key="col.name" :props="props">
                                {{ col.label }}
                            </q-th>
                        </q-tr>
                    </template>

                    <template v-slot:body="props">
                        <q-tr :props="props">
                            <q-td auto-width>
                                <q-btn size="sm" color="accent" round dense @click="props.expand = !props.expand"
                                    :icon="props.expand ? 'remove' : 'add'" />
                            </q-td>
                            <q-td v-for="col in props.cols" :key="col.name" :props="props">
                                {{ col.value }}
                            </q-td>
                        </q-tr>
                        <q-tr v-show="props.expand" :props="props">
                            <q-td colspan="100%">
                                <div class="text-caption">MSA prediction</div>
                                <q-skeleton :type="text" />
                                <div class="text-caption">Folding</div>
                                <q-skeleton :type="text" />
                                <div class="text-caption">Pocket prediction</div>
                                <q-skeleton :type="text" />
                                <div class="text-caption">Docking</div>
                                <q-skeleton :type="text" />
                            </q-td>
                        </q-tr>
                    </template>


                </q-table>
            </div>
        </q-page>
    </q-page>
</template>

<script>
import { useDrugDiscoveryStore } from 'src/features/drug_discovery/storage';
import { QSpinner } from 'quasar';
import { useRoute } from 'vue-router';

export default {
    name: 'RunningJobs',
    components: {
    },
    data() {
        return {
            loading: true,
            running_jobs: [
                {
                    id: 12312345,
                    target: '7nb4',
                    ligand: 'donepizil'
                },
                {
                    id: 1231234,
                    target: '7nb4',
                    ligand: 'ligand2'
                },
                {
                    id: 12312342,
                    target: '7nb4',
                    ligand: 'ligand3'
                }
            ],
            queue_jobs: [
                {
                    id: 1231234,
                    target: '7nb4',
                    ligand: 'donepizil'
                },
                {
                    id: 1231234,
                    target: '7nb4',
                    ligand: 'ligand2'
                },
                {
                    id: 1231234,
                    target: '7nb4',
                    ligand: 'ligand3'
                }
            ],
            columns: [
                { name: 'id', align: 'left', label: 'Id', field: 'id', sortable: true },
                { name: 'target', align: 'left', label: 'Target', field: 'target' },
                { name: 'carbs', align: 'left', label: 'Ligand', field: 'ligand' },
            ]
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
