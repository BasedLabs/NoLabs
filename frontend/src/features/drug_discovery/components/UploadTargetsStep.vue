<template>
  <!-- Main Page Content -->
  <q-page>
    <div v-if="loading">
      <q-spinner color="primary"/>
    </div>
    <div v-else>
      <TargetsList
          :experimentId="this.experimentId"
          :originalTargets="this.targets"
      />
      <div class="q-pt-xl"></div>
      <div class="q-pt-xl"></div>
      <LigandsList
        v-if="ligands"
        :experimentId="this.experimentId"
        :originalLigands="this.ligands"
      />
      <q-expansion-item
          icon="info"
          label="Step 1: Upload target files for which you discover the drugs"
          class="q-mt-md"
          expand-icon="eye"
          expanded-icon="minimize"
          :dense="true"
      >
        <q-card-section class="q-pa-md">
          <p>
            Upload the <strong>.fasta</strong> files of the proteins to which you want to design drugs. If you don't
            have protein files prepared,
            you can visit <a style="color: #31ccec" href="https://www.rcsb.org/"
                             target="_blank">https://www.rcsb.org/</a> to download .fasta files of your interest.
          </p>
          <p class="q-mt-md">
            <strong>Warning:</strong> If you have a light infrastructure set in your <code>nolabs/infrastructure/settings.ini</code>
            file, then sequences longer than 400 amino acids will not be processed.
          </p>
          <p class="q-mt-md">
            Binding pockets can be set via UI when a user clicks on the target, or they will be predicted automatically.
          </p>
          <div class="q-pa-md q-gutter-sm row justify-center">
            <q-img :src="'/drugDiscovery/Protein_grey_1.png'" class="col" style="max-width: 200px"/>
            <q-img :src="'/drugDiscovery/Protein_grey_2.png'" class="col" style="max-width: 200px"/>
          </div>
        </q-card-section>
      </q-expansion-item>
    </div>
  </q-page>
</template>

<script lang="ts">
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {Notify, QSpinnerOrbit} from "quasar";
import {useRoute} from "vue-router";
import TargetsList from "src/features/drug_discovery/components/targets/TargetsList.vue";
import {defineComponent} from "vue";
import {LigandMetaData, TargetMetaData} from "../../../api/client";
import LigandsList from "./ligands/LigandsList.vue";
import {ExtendedLigandMetaData, ExtendedTargetMetaData} from "./targets/types";

export default defineComponent({
  name: "ExperimentSetup",
  components: {
    TargetsList,
    LigandsList,
  },
  data() {
    return {
      loading: true,
      experimentId: null as string | null
    };
  },
  computed: {
    targets(): ExtendedTargetMetaData[] {
      const store = useDrugDiscoveryStore();
      return store.targets;
    },
    ligands(): ExtendedLigandMetaData[] {
      const store = useDrugDiscoveryStore();
      return store.loneLigands;
    },
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
      await store.fetchLigandsForExperiment(this.experimentId);
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
