<template>
  <ExperimentsView :store="store" path-to-experiment-page="Folding"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import useFoldingStore from "src/features/aminoAcid/folding/storage";
import ExperimentsView from "src/components/ExperimentsView.vue";


export default defineComponent({
  name: 'FoldingExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useFoldingStore()
    }
  }
})
</script>
