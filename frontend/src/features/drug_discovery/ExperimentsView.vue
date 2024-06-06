<template>
  <ExperimentsView :store="store" path-to-experiment-page="Drug discovery"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import { useWorkflowStore } from 'src/features/drug_discovery/components/workflow/storage';


export default defineComponent({
  name: 'DrugTargetExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useWorkflowStore()
    }
  }
})
</script>
