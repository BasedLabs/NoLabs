<template>
    <div class="q-pa-md">
        <q-table flat bordered grid title="Experiments" v-model:selected="selected" :rows="rows" :columns="columns"
            :loading="loading" row-key="id" no-data-label="No experiments found. Add a new experiment." :filter="filter"
            hide-header>
            <template v-slot:top-right>
                <q-btn @click="addExperimentClick" class="q-mx-md" color="positive" icon-right="add" label="Add experiment"
                    no-caps />
                <q-input borderless dense debounce="300" v-model="filter" placeholder="Search">
                    <template v-slot:append>
                        <q-icon name="search" />
                    </template>
                </q-input>
            </template>
            <template v-slot:item="props">
                <router-link style="text-decoration: none;"
                    :to="{ path: `/labs/protein-design/experiment/${props.row.id}` }">
                    <div class="q-pa-xs col-xs-12 col-sm-6 col-md-4 col-lg-3 grid-style-transition"
                        :style="props.selected ? 'transform: scale(0.95);' : ''">
                        <q-card bordered flat :class="props.selected ? ($q.dark.isActive ? 'bg-grey-9' : 'bg-grey-2') : ''">
                            <q-list dense>
                                <q-item v-for="col in props.cols.filter(col => col.name !== 'desc')" :key="col.name">
                                    <q-item-section>
                                        <q-item-label>{{ col.label }}</q-item-label>
                                    </q-item-section>
                                    <q-item-section side>
                                        <q-item-label caption>{{ col.value }}</q-item-label>
                                    </q-item-section>
                                </q-item>
                            </q-list>
                        </q-card>
                    </div>
                </router-link>
            </template>
        </q-table>
    </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import useProteinDesignStore from './storage';
import { ExperimentListItem } from './types';


export default defineComponent({
    name: 'ExperimentsView',
    data() {
        const store = useProteinDesignStore();
        const columns = [
            {
                name: 'name',
                label: 'Experiment Name',
                field: row => row.name,
                style: 'font-size: 3em'
            }, {
                name: 'id',
                label: 'Experiment Id',
                field: row => row.id,
                style: 'font-size: 3em'
            }];
        let experiments: ExperimentListItem[] = [];
        return {
            store,
            columns,
            rows: experiments,
            filter: '',
            selected: [],
            loading: false,
            initialPagination: {
                sortBy: 'desc',
                descending: false,
                page: 1,
                rowsPerPage: 10
                // rowsNumber: xx if getting data from a server
            },
        }
    },
    methods: {
        addExperimentClick() {
            this.store.addExperiment();
        },
    },
    async mounted() {
        this.loading = true;
        await this.store.getExperiments();
        this.rows = this.store.experiments;
        this.loading = false;
    }
})
</script>