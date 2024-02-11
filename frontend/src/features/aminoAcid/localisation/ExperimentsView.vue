<template>
  <ExperimentsView :store="store" path-to-experiment-page="Localisation"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import useLocalisationStore from "src/features/aminoAcid/localisation/storage";


export default defineComponent({
  name: 'LocalisationExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useLocalisationStore()
    }
  }
})
</script>
