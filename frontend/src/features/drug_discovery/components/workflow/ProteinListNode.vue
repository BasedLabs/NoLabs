<template>
    <q-card>
      <TargetsList v-if="experimentId" :experiment-id="this.experimentId" />
    </q-card>
</template>

<script lang="ts">
import {defineComponent} from "vue";
import TargetsList from "../targets/TargetsList.vue";
import {useDrugDiscoveryStore} from "../../storage";
import {useRoute} from "vue-router";
import {Notify, QSpinnerOrbit} from "quasar";

export default defineComponent({
  name: "LigandsListNode",
  components: {
    TargetsList,
  },
  data() {
    return {
      loading: true,
      experimentId: null as string | null,
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    const route = useRoute();

    this.experimentId = route.params.experimentId as string;
    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Fetching targets for experiment`
    });
    try {
      await store.fetchTargetsForExperiment(this.experimentId);
    } catch (e) {
      Notify.create({
        type: "negative",
        closeBtn: 'Close',
        message: 'Error fetching targets or ligands: ' + e
      });
    } finally {
      this.loading = false;
      this.$q.loading.hide();
    }
  },
})
</script>

<style scoped>
.special-node {
  border: 2px solid #ff6600;
  background-color: black;
  padding: 8px;
  border-radius: 4px;
  height: 400px;
  width: 400px
}

.special-node.input {
  border-color: #00bcd4; /* Change border color for input nodes */
}

.special-node.output {
  border-color: #4caf50; /* Change border color for output nodes */
}

.content {
  height: 100%;
  width: 100%;
}
</style>
