<template>
  <ExperimentsView :store="store" path-to-experiment-page="Small molecules design"/>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import ExperimentsView from "src/components/ExperimentsView.vue";
import useSmallMoleculesDesignStore from "./storage";


export default defineComponent({
  name: 'SmallMoleculeDesignExperimentsView',
  components: {ExperimentsView},
  data() {
    return {
      store: useSmallMoleculesDesignStore()
    }
  }
})
</script>
