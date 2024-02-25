<template>
  <div class="umol-result">
    <q-page padding>
      <PdbViewer v-if="predictedPdb && predictedSDF" :pdb-file="predictedPdb" :sdf-file="predictedSDF" />

      <div v-if="plddtArray && plddtArray.length" class="q-mt-md">
        <h3>PLDDT Scores:</h3>
        <div class="plddt-scores q-gutter-sm q-mb-md">
          <q-chip v-for="(score, index) in plddtArray" :key="index" outline dense>
            {{ score }}
          </q-chip>
        </div>
      </div>

      <div v-if="pocketIds && pocketIds.length" class="q-mt-md">
        <h3>Pocket IDs:</h3>
        <div class="pocket-ids q-gutter-sm">
          <q-chip v-for="id in pocketIds" :key="id" outline color="primary" dense>
            {{ id }}
          </q-chip>
        </div>
      </div>
    </q-page>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QPage, QChip } from 'quasar';
import PdbViewer from 'src/components/PdbViewer.vue';
import {useDrugDiscoveryStore} from 'src/features/drug_discovery/storage';

export default defineComponent({
  name: 'UmolResult',
  components: {
    PdbViewer,
    QPage,
    QChip,
  },
  props: {
    experimentId: String,
    targetId: String,
    ligandId: String,
    jobId: String,
  },
  data() {
    return {
      predictedPdb: null as File | null,
      predictedSDF: null as File | null,
      plddtArray: [] as number[],
      pocketIds: [] as number[],
    };
  },
  mounted() {
    this.loadUmolResults();
  },
  methods: {
    async loadUmolResults() {
      const store = useDrugDiscoveryStore();
      const resultData = await store.getUmolDockingJobResultData(this.experimentId!, this.targetId!, this.ligandId!, this.jobId!);
      if (resultData) {
        this.predictedPdb = new File([new Blob([resultData.predicted_pdb])], "predicted.pdb");
        this.predictedSDF = new File([new Blob([resultData.predicted_sdf])], "predicted.sdf");
        this.plddtArray = resultData.plddt_array;
      }
      const pocketData = await store.getJobPocketIds(this.experimentId!, this.targetId!, this.ligandId!, this.jobId!);
      if (pocketData && pocketData.pocket_ids) {
        this.pocketIds = pocketData.pocket_ids;
      }
    },
  },
});
</script>

<style scoped>
.plddt-scores, .pocket-ids {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
}
</style>
