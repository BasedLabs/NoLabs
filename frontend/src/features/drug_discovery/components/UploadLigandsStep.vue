<template>
  <!-- Main Page Content -->
  <q-page class="q-ma-md">
    <div v-if="loading">
      <q-spinner color="primary" />
    </div>
    <div v-else>
      <TargetsLigandsList
        :experimentId="this.experimentId"
        :originalTargets="this.targets"
      />
    </div>
  </q-page>
</template>

<script>
import { useDrugDiscoveryStore } from "src/features/drug_discovery/storage";
import { QSpinner } from "quasar";
import { useRoute } from "vue-router";
import TargetsLigandsList from "src/features/drug_discovery/components/targets/TargetsLigandsList.vue";

export default {
  name: "ExperimentSetup",
  components: {
    TargetsLigandsList,
  },
  data() {
    return {
      loading: true,
      targets: null,
    };
  },
  mounted() {
    const store = useDrugDiscoveryStore();
    const route = useRoute();
    this.experimentId = route.params.experimentId;
    store
      .fetchTargetsForExperiment(this.experimentId)
      .then(() => {
        this.targets = store.targets;
        this.loading = false;
      })
      .catch((error) => {
        console.error("Error fetching targets:", error);
        this.loading = false;
      });
  },
};
</script>
