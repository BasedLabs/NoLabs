<template>
  <q-card>
    <TargetsList v-if="experimentId" :experiment-id="experimentId" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import TargetsList from '../targets/TargetsList.vue';
import { useDrugDiscoveryStore } from '../../storage';
import { useRoute } from 'vue-router';
import { Notify, QSpinnerOrbit } from 'quasar';

export default defineComponent({
  name: 'ProteinListNodeContent',
  components: { TargetsList },
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
      message: `Fetching targets for experiment`,
    });
    try {
      await store.fetchTargetsForExperiment(this.experimentId);
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
