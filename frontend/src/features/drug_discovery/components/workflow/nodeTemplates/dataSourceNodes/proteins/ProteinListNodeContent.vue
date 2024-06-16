<template>
    <q-card>
      <ProteinList v-if="experimentId" :experiment-id="experimentId" :nodeId="nodeId" />
    </q-card>
  </template>
  
  <script lang="ts">
  import { defineComponent } from 'vue';
  import { Notify, QSpinnerOrbit } from 'quasar';
  import ProteinList from './ProteinList.vue';
  
  export default defineComponent({
    name: 'ProteinListNodeContent',
    components: { ProteinList },
    props: {
      experimentId: {
            type: String,
            required: true,
        },
        nodeId: {
            type: String,
            required: true,
        }
    },
    data() {
      return {
        loading: true,
      };
    },
    async mounted() {
    
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Fetching targets for experiment`,
      });
      try {
        //await store.fetchTargetsForExperiment(this.experimentId);
      } catch (e) {
        Notify.create({
          type: 'negative',
          closeBtn: 'Close',
          message: 'Error fetching targets: ' + e,
        });
      } finally {
        this.loading = false;
        this.$q.loading.hide();
      }
    },
  });
  </script>
  
  <style scoped>
  </style>
  