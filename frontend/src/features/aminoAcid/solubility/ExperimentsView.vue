<template>
  <ExperimentsView :store="store" path-to-experiment-page="Solubility"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import useSolubilityStore from "src/features/aminoAcid/solubility/storage";


export default defineComponent({
  name: 'SolubilityExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useSolubilityStore()
    }
  }
})
</script>
