<script lang="ts">
import { defineComponent, PropType } from "vue";

type Experiment = {
  id: string;
  name: string;
};

export { Experiment };

export default defineComponent({
  name: "ExperimentsList",
  props: {
    experiments: [] as PropType<Array<Experiment>>,
    onExperimentClick: ((experiment: Experiment) => {}) as PropType<any>,
  },
  data() {
    activeExperiment: null as Experiment | null;
  },
  methods: {
    experimentClick(experiment: Experiment) {
      this.activeExperiment = experiment;
      this.onExperimentClick(experiment);
    },
  },
});
</script>

<template>
  ASDASDASD
  <q-list padding class="menu-list">
    <q-item
      v-on:click="this.experimentClick(experiment)"
      clickable
      :key="experiment.id"
      :active="experiment === this.activeExperiment"
      v-ripple
      v-for="experiment in this.experiments"
    >
      <q-item-section>
        {{ experiment.name }}
      </q-item-section>
    </q-item>
  </q-list>
</template>
