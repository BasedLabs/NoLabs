<template>
  <ExperimentsView :store="store" path-to-experiment-page="DrugDiscoveryExperiment"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";


export default defineComponent({
  name: 'DrugTargetExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useDrugDiscoveryStore()
    }
  }
})
</script>
