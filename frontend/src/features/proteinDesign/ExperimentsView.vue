<template>
  <ExperimentsView :store="store" path-to-experiment-page="Protein design"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import useProteinDesignStore from "src/features/proteinDesign/storage";


export default defineComponent({
  name: 'SolubilityExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useProteinDesignStore()
    }
  }
})
</script>
