<template>
    <!-- Main Page Content -->
    <q-page class="q-ma-md">
      <div v-if="loading">
        <q-spinner color="primary" />
      </div>
      <div v-else>
        <TargetsList
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
import TargetsList from "src/features/drug_discovery/components/targets/TargetsList.vue";

export default {
  name: "ExperimentSetup",
  components: {
    TargetsList,
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
