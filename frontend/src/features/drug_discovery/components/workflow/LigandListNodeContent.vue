<template>
  <q-card>
    <LigandsList v-if="experimentId" :experiment-id="experimentId" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import LigandsList from '../ligands/LigandsList.vue';
import { useDrugDiscoveryStore } from '../../storage';
import { useRoute } from 'vue-router';
import { Notify, QSpinnerOrbit } from 'quasar';

export default defineComponent({
  name: 'LigandsNodeContent',
  components: { LigandsList },
  props: {
    experimentId: String,
  },
  data() {
    return {
      loading: true,
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Fetching ligands for experiment`,
    });
    try {
      await store.fetchLigandsForExperiment(this.experimentId);
    } catch (e) {
      Notify.create({
        type: 'negative',
        closeBtn: 'Close',
        message: 'Error fetching ligands: ' + e,
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
