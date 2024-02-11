<template>
  <!-- Main Page Content -->
  <q-page class="q-ma-md">
    <div v-if="loading">
      <q-spinner color="primary"/>
    </div>
    <div v-else>
      <TargetsLigandsList
          :experimentId="this.experimentId"
          :originalTargets="this.targets"
      />
      <q-expansion-item
          icon="info"
          label="Step 2: Upload ligand (drug) .sdf files"
          class="q-mt-md"
          expand-icon="eye"
          expanded-icon="minimize"
          :dense="true"
      >
        <q-card-section class="q-pa-md">
          <p>
            Upload the <strong>.sdf</strong> files of ligands (drugs) to targets of interest. If you don't have ligand
            files prepared,
            you can visit <a style="color: #31ccec" href="https://www.rcsb.org/downloads/ligands" target="_blank">https://www.rcsb.org/downloads/ligands</a>
            to download .sdf files of your interest.
          </p>
          <div class="q-pa-md q-gutter-sm row justify-center">
            <q-img :src="'/drugDiscovery/Upload_ligands.png'" class="col" style="max-width: 400px"/>
          </div>
        </q-card-section>
      </q-expansion-item>
    </div>
  </q-page>
</template>

<script>
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {Notify, QSpinner, QSpinnerOrbit} from "quasar";
import {useRoute} from "vue-router";
import TargetsLigandsList from "src/features/drug_discovery/components/targets/TargetsLigandsList.vue";

export default {
  name: "ExperimentSetup",
  components: {
    TargetsLigandsList,
  },
  data() {
    return {
      loading: true,
      targets: null
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    const route = useRoute();
    this.experimentId = route.params.experimentId;
    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Fetching targets for experiment`
    });
    try {
      await store.fetchTargetsForExperiment(this.experimentId);
      this.targets = store.targets;
    } catch (e) {
      Notify.create({
        type: "negative",
        closeBtn: 'Close',
        message: 'Error fetching targets: ' + e
      });
    } finally {
      this.loading = false;
      this.$q.loading.hide();
    }
  },
};
</script>
