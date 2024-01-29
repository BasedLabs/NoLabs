<template>
  <ExperimentsView :store="store" path-to-experiment-page="ConformationsExperiment"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import useConformationsStore from "src/features/conformations/storage";


export default defineComponent({
  name: 'SolubilityExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useConformationsStore()
    }
  }
})
</script>
