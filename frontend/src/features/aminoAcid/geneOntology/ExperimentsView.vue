<template>
  <ExperimentsView :store="store" path-to-experiment-page="Gene ontology"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import useGeneOntologyStore from "src/features/aminoAcid/geneOntology/storage";


export default defineComponent({
  name: 'GeneOntologyExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useGeneOntologyStore()
    }
  }
})
</script>
