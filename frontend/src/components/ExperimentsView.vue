<template>
  <div class="q-pa-md">
    <q-table flat grid bordered title="Experiments" v-model:selected="selected" :rows="experiments" :columns="columns"
             :loading="loading" row-key="id" no-data-label="No experiments found. Add a new experiment."
             :filter="filter"
             hide-header>
      <template v-slot:top-right>
        <q-btn @click="addExperimentClick()" outline class="q-mx-md" color="info" icon-right="add"
               label="Add experiment" no-caps/>
        <q-input borderless dense debounce="300" v-model="filter" placeholder="Search">
          <template v-slot:append>
            <q-icon name="search"/>
          </template>
        </q-input>
      </template>
      <template v-slot:item="props">
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
            <q-card-actions align="right">
              <router-link style="text-decoration: none;"
                           :to="{ name: 'Workflow view', params: { experimentId: props.row.id } }">
                <q-btn outline color="info" label="Open"/>
              </router-link>
              <q-btn outline class="q-mx-md" color="negative" label="Remove"
                     @click="removeExperimentClick(props.row.id)"/>
            </q-card-actions>
          </q-card>
        </div>
      </template>
    </q-table>
  </div>
</template>

<script lang="ts">
import {defineComponent, PropType} from 'vue';
import {ExperimentListItem} from "src/components/types";
import useExperimentsStore from "./experimentsStorage";


export default defineComponent({
  name: 'ExperimentsView',
  data() {
    return {
      experiments: [] as ExperimentListItem[],
      filter: '',
      selected: [],
      loading: false,
      initialPagination: {
        sortBy: 'desc',
        descending: false,
        page: 1,
        rowsPerPage: 10
      }
    }
  },
  store: useExperimentsStore(),
  computed: {
    columns() {
      return [
        {
          name: 'name',
          label: 'Experiment Name',
          field: (row: {name: string, id:string}) => row.name,
          style: 'font-size: 3em'
        }, {
          name: 'id',
          label: 'Experiment Id',
          field: (row: {name: string, id:string}) => row.id,
          style: 'font-size: 3em'
        }
      ]
    }
  },
  methods: {
    async addExperimentClick() {
      const response = await this.$options.store.createExperiment();
      this.experiments.push({
        id: response!.id,
        name: response!.name
      })
    },
    async removeExperimentClick(experimentId: string) {
      await this.$options.store.deleteExperiment(experimentId);
      this.experiments = this.experiments.filter(x => x.id !== experimentId);
    }
  },
  async mounted() {
    this.loading = true;
    const response = await this.$options.store.getExperiments();
    if(response.length > 0){
      this.experiments = response;
    }
    this.loading = false;
  }
})
</script>
